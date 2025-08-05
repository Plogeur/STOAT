#include <stdio.h>  // printf
#include <algorithm>  // fill
#include <cassert>
#include <numeric>  // accumulate

#include "third/eigen/Eigen/Cholesky"
#include "third/gsl/include/gsl/gsl_cdf.h"
#include "third/eigen/Eigen/Core"

void Vector::Dimension(int n) { data.resize(n); }
void Vector::Dimension(int n, double val) {
  data.resize(n);
  Fill(val);
}

void Vector::Fill(double val) { std::fill(data.begin(), data.end(), val); }
double Vector::Sum() const {
  return std::accumulate(data.begin(), data.end(), 0.0);
}
double Vector::Average() const {
  if (data.empty()) return 0.0;
  return Sum() / data.size();
}
double Vector::Min() const {
  return *std::min_element(data.begin(), data.end());
}
double Vector::Max() const {
  return *std::max_element(data.begin(), data.end());
}

Matrix::Matrix(int nr, int nc) {
  rows = nr;
  cols = nc;
  data.resize(nr * nc);
}

Matrix::Matrix(const Matrix& m) {
  rows = m.rows;
  cols = m.cols;
  data = m.data;
  colLabel = m.colLabel;

  static int i = 0;
  ++i;
  printf("Matrix() called %d times\n", i);
}

Matrix& Matrix::operator=(const Matrix& m) {
  rows = m.rows;
  cols = m.cols;
  data = m.data;
  colLabel = m.colLabel;

  return *this;
}

void Matrix::Dimension(int nr, int nc) {
  if (nr == rows && nc == cols) {
    return;
  }

  // todo: make this run faster
  std::vector<double> newData(nr * nc);
  for (int i = 0; i < nr && i < rows; ++i) {
    for (int j = 0; j < nc && j < cols; ++j) {
      newData[i + j * nr] = data[i + j * rows];
    }
  }

  rows = nr;
  cols = nc;
  std::swap(data, newData);
  colLabel.resize(nc);
}

/**
 * Set all matrix elements to @param val
 */
void Matrix::Dimension(int nr, int nc, double val) {
  DimensionQuick(nr, nc);
  Fill(val);
}

void Matrix::DimensionQuick(int nr, int nc) {
  assert(nr >= 0 && nc >= 0);
  rows = nr;
  cols = nc;
  data.resize(nr * nc);
  colLabel.resize(nc);
}

void Matrix::Reserve(int nr, int nc) {
  data.reserve(nr * nc);
  colLabel.reserve(nc);
}

double Matrix::Min() const {
  return *std::min_element(data.begin(), data.end());
}

double Matrix::Max() const {
  return *std::max_element(data.begin(), data.end());
}

void Matrix::Product(const Matrix& in1, const Matrix& in2) {
  DECLARE_EIGEN_CONST_MATRIX(in1, in1_e);
  DECLARE_EIGEN_CONST_MATRIX(in2, in2_e);
  DimensionQuick(in1_e.rows(), in2_e.cols());
  Eigen::Map<Eigen::MatrixXd> ret(data.data(), rows, cols);
  ret = in1_e * in2_e;
}

void Matrix::Transpose(const Matrix& old) {
  data.resize(old.data.size());
  rows = old.cols;
  cols = old.rows;
  DECLARE_EIGEN_CONST_MATRIX(old, old_e);
  DECLARE_EIGEN_MATRIX((*this), new_e);
  new_e = old_e.transpose();
}

Matrix& Matrix::Multiply(double s) {
  for (std::vector<double>::iterator iter = data.begin(); iter != data.end();
       ++iter) {
    *iter *= s;
  }
  return *this;
}

int Matrix::RemoveByRowIndex(const std::vector<int>& rowIndexToRemove) {
  int idx = 0;
  assert(*std::min_element(rowIndexToRemove.begin(), rowIndexToRemove.end()) >=
         0);
  assert(*std::max_element(rowIndexToRemove.begin(), rowIndexToRemove.end()) <
         rows);
  std::set<int> idxSet(rowIndexToRemove.begin(), rowIndexToRemove.end());
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      if (idxSet.count(i)) {
        continue;
      }
      data[idx++] = (*this)(i, j);
    }
  }
  rows -= idxSet.size();
  data.resize(rows * cols);
  return idxSet.size();
}

Matrix& Matrix::StackRight(const Matrix& m) {
  assert(rows = m.rows);
  data.insert(data.end(), m.data.begin(), m.data.end());
  cols += m.cols;
  colLabel.insert(colLabel.end(), m.colLabel.begin(), m.colLabel.end());
  return *this;
}

bool LinearRegression::FitLinearModel(const Matrix& X, const Vector& y) {
  // Matrix Xt;
  // Xt.Transpose(X);

  // Matrix XtX;
  // XtX.Product(Xt, X);
  // if (!this->chol.TryDecompose(XtX)) return false;
  // chol.Decompose(XtX);
  // chol.Invert();
  // this->XtXinv = chol.inv;
  XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e)
                 .llt()
                 .solve(Eigen::MatrixXd::Identity(X_e.cols(), X_e.cols()));

  // Vector tmp = y;
  // tmp.Product(Xt, y);
  // this->B.Product(this->XtXinv, tmp);  // beta = (XtX)^{-1} Xt Y
  B.Dimension(X.cols, 1);
  DECLARE_EIGEN_VECTOR(B, B_e);
  DECLARE_EIGEN_CONST_VECTOR(y, y_e);
  B_e = XtXinv_e * X_e.transpose() * y_e;

  // this->predict.Product(X, this->B);
  // this->residuals = y;
  // this->residuals.Subtract(this->predict);
  this->predict.Dimension(X.rows, 1);
  this->residuals.Dimension(X.rows, 1);
  DECLARE_EIGEN_VECTOR(this->predict, predict_e);
  DECLARE_EIGEN_VECTOR(this->residuals, resid_e);
  predict_e = X_e * B_e;
  resid_e = y_e - predict_e;

  // this->sigma2 = 0.0;
  // for (int i = 0; i < this->residuals.Length(); i++) {
  //   sigma2 += (this->residuals[i]) * (this->residuals[i]);
  // }
  // sigma2 /= y.Length();  // MLE estimates of sigma2
  this->sigma2 = resid_e.squaredNorm() / y_e.size();

  // this->covB = this->XtXinv;
  // this->covB.Multiply(sigma2);
  this->covB.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(this->covB, covB_e);
  covB_e = XtXinv_e * sigma2;

  return true;
};

Vector& LinearRegression::GetAsyPvalue() {
  int numCov = B.Length();
  pValue.Dimension(B.Length());
  for (int i = 0; i < numCov; i++) {
    double Zstat = B[i] / sqrt(covB(i, i));
    // pValue[i] = ndist(Zstat);
    // if (pValue[i] >= 0.5){
    //      pValue[i] = 2*(1-pValue[i]);
    // } else pValue[i] = 2*pValue[i];
    Zstat *= Zstat;
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return (pValue);
}

bool LinearRegression::calculateResidualMatrix(Matrix& X, Matrix* out) {
  if (!calculateHatMatrix(X, out)) return false;
  Matrix& m = *out;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      if (i == j) {
        m(i, j) = 1.0 - m(i, j);
      } else {
        m(i, j) = -m(i, j);
      }
    }
  }
  return true;
}

bool LinearRegression::calculateHatMatrix(Matrix& X, Matrix* out) {
  // Matrix Xt;
  // Xt.Transpose(X);

  // Matrix XtX;
  // XtX.Product(Xt, X);
  // if (!this->chol.TryDecompose(XtX)) return false;
  // chol.Decompose(XtX);
  // chol.Invert();
  // this->XtXinv = chol.inv;
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  this->XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(this->XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e)
                 .llt()
                 .solve(Eigen::MatrixXd::Identity(X_e.cols(), X_e.cols()));

  // Matrix tmp;
  // tmp.Product(XtXinv, Xt);
  // (*out).Product(X, tmp);

  DECLARE_EIGEN_MATRIX((*out), out_e);
  out_e = X_e * XtXinv_e * X_e.transpose();
  return true;
}