#include <stdio.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <set>
#include <vector>

#include "third/eigen/Eigen/Cholesky"
#include "third/eigen/Eigen/Core"
#include "third/gsl/include/gsl/gsl_cdf.h"

#define DECLARE_EIGEN_VECTOR(v, v_e) Eigen::Map<Eigen::VectorXd> v_e((v).data.data(), (v).data.size())
#define DECLARE_EIGEN_CONST_VECTOR(v, v_e) Eigen::Map<const Eigen::VectorXd> v_e((v).data.data(), (v).data.size())
#define DECLARE_EIGEN_MATRIX(m, m_e) Eigen::Map<Eigen::MatrixXd> m_e((m).data.data(), (m).rows, (m).cols)
#define DECLARE_EIGEN_CONST_MATRIX(m, m_e) Eigen::Map<const Eigen::MatrixXd> m_e((m).data.data(), (m).rows, (m).cols)

// ======================= Vector Class =========================
class Vector {
 public:
  std::vector<double> data;

  Vector() {}
  Vector(int n) { Dimension(n); }
  double& operator[](int i) { return data[i]; }
  double operator[](int i) const { return data[i]; }
  int Length() const { return data.size(); }

  void Dimension(int n);
  void Dimension(int n, double val);
  void Fill(double val);
  double Sum() const;
  double Average() const;
  double Min() const;
  double Max() const;
};

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
  return data.empty() ? 0.0 : Sum() / data.size();
}
double Vector::Min() const {
  return *std::min_element(data.begin(), data.end());
}
double Vector::Max() const {
  return *std::max_element(data.begin(), data.end());
}

// ======================= Matrix Class =========================
class Matrix {
 public:
  int rows, cols;
  std::vector<double> data;
  std::vector<std::string> colLabel;

  Matrix(int nr, int nc);
  Matrix(const Matrix& m);
  Matrix& operator=(const Matrix& m);

  double& operator()(int r, int c) { return data[r + c * rows]; }
  double operator()(int r, int c) const { return data[r + c * rows]; }

  void Dimension(int nr, int nc);
  void Dimension(int nr, int nc, double val);
  void DimensionQuick(int nr, int nc);
  void Reserve(int nr, int nc);
  void Fill(double val);
  double Min() const;
  double Max() const;

  void Product(const Matrix& in1, const Matrix& in2);
  void Transpose(const Matrix& old);
  Matrix& Multiply(double s);
  int RemoveByRowIndex(const std::vector<int>& rowIndexToRemove);
  Matrix& StackRight(const Matrix& m);
};

Matrix::Matrix(int nr, int nc) : rows(nr), cols(nc), data(nr * nc) {}

Matrix::Matrix(const Matrix& m)
    : rows(m.rows), cols(m.cols), data(m.data), colLabel(m.colLabel) {}

Matrix& Matrix::operator=(const Matrix& m) {
  rows = m.rows;
  cols = m.cols;
  data = m.data;
  colLabel = m.colLabel;
  return *this;
}

void Matrix::Dimension(int nr, int nc) {
  if (nr == rows && nc == cols) return;
  std::vector<double> newData(nr * nc);
  for (int i = 0; i < nr && i < rows; ++i)
    for (int j = 0; j < nc && j < cols; ++j)
      newData[i + j * nr] = data[i + j * rows];
  rows = nr;
  cols = nc;
  std::swap(data, newData);
  colLabel.resize(nc);
}

void Matrix::Dimension(int nr, int nc, double val) {
  DimensionQuick(nr, nc);
  Fill(val);
}

void Matrix::DimensionQuick(int nr, int nc) {
  rows = nr;
  cols = nc;
  data.resize(nr * nc);
  colLabel.resize(nc);
}

void Matrix::Reserve(int nr, int nc) {
  data.reserve(nr * nc);
  colLabel.reserve(nc);
}

void Matrix::Fill(double val) {
  std::fill(data.begin(), data.end(), val);
}

double Matrix::Min() const {
  return *std::min_element(data.begin(), data.end());
}

double Matrix::Max() const {
  return *std::max_element(data.begin(), data.end());
}

void Matrix::Product(const Matrix& in1, const Matrix& in2) {
  DECLARE_EIGEN_CONST_MATRIX(in1, e1);
  DECLARE_EIGEN_CONST_MATRIX(in2, e2);
  DimensionQuick(in1.rows, in2.cols);
  DECLARE_EIGEN_MATRIX((*this), out);
  out = e1 * e2;
}

void Matrix::Transpose(const Matrix& old) {
  DimensionQuick(old.cols, old.rows);
  DECLARE_EIGEN_CONST_MATRIX(old, eOld);
  DECLARE_EIGEN_MATRIX((*this), eNew);
  eNew = eOld.transpose();
}

Matrix& Matrix::Multiply(double s) {
  for (double& val : data) val *= s;
  return *this;
}

int Matrix::RemoveByRowIndex(const std::vector<int>& rowIndexToRemove) {
  std::set<int> idxSet(rowIndexToRemove.begin(), rowIndexToRemove.end());
  int idx = 0;
  for (int j = 0; j < cols; ++j)
    for (int i = 0; i < rows; ++i)
      if (!idxSet.count(i)) data[idx++] = (*this)(i, j);
  rows -= idxSet.size();
  data.resize(rows * cols);
  return idxSet.size();
}

Matrix& Matrix::StackRight(const Matrix& m) {
  assert(rows == m.rows);
  data.insert(data.end(), m.data.begin(), m.data.end());
  cols += m.cols;
  colLabel.insert(colLabel.end(), m.colLabel.begin(), m.colLabel.end());
  return *this;
}

// ==================== LinearRegression Class ===================
class LinearRegression {
 public:
  Matrix XtXinv, B, covB;
  Vector predict, residuals, pValue;
  double sigma2;

  bool FitLinearModel(const Matrix& X, const Vector& y);
  Vector& GetAsyPvalue();
  bool calculateHatMatrix(Matrix& X, Matrix* out);
  bool calculateResidualMatrix(Matrix& X, Matrix* out);
};

bool LinearRegression::FitLinearModel(const Matrix& X, const Vector& y) {
  XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e).llt().solve(Eigen::MatrixXd::Identity(X.cols, X.cols));

  B.Dimension(X.cols, 1);
  DECLARE_EIGEN_VECTOR(B, B_e);
  DECLARE_EIGEN_CONST_VECTOR(y, y_e);
  B_e = XtXinv_e * X_e.transpose() * y_e;

  predict.Dimension(X.rows, 1);
  residuals.Dimension(X.rows, 1);
  DECLARE_EIGEN_VECTOR(predict, pred_e);
  DECLARE_EIGEN_VECTOR(residuals, resid_e);
  pred_e = X_e * B_e;
  resid_e = y_e - pred_e;

  sigma2 = resid_e.squaredNorm() / y.Length();

  covB.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(covB, covB_e);
  covB_e = XtXinv_e * sigma2;

  return true;
}

Vector& LinearRegression::GetAsyPvalue() {
  int numCov = B.Length();
  pValue.Dimension(numCov);
  for (int i = 0; i < numCov; ++i) {
    double Zstat = B[i] / sqrt(covB(i, i));
    Zstat *= Zstat;
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return pValue;
}

bool LinearRegression::calculateHatMatrix(Matrix& X, Matrix* out) {
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e).llt().solve(Eigen::MatrixXd::Identity(X.cols, X.cols));
  DECLARE_EIGEN_MATRIX((*out), out_e);
  out_e = X_e * XtXinv_e * X_e.transpose();
  return true;
}

bool LinearRegression::calculateResidualMatrix(Matrix& X, Matrix* out) {
  if (!calculateHatMatrix(X, out)) return false;
  for (int i = 0; i < out->rows; ++i)
    for (int j = 0; j < out->cols; ++j)
      (*out)(i, j) = (i == j) ? 1.0 - (*out)(i, j) : -(*out)(i, j);
  return true;
}

// ======================= Test Main ===========================
int main() {
  // Simple test with 3 data points and 2 variables (X0 = 1 for intercept)
  Matrix X(3, 2);
  X(0, 0) = 1; X(0, 1) = 1;
  X(1, 0) = 1; X(1, 1) = 2;
  X(2, 0) = 1; X(2, 1) = 3;

  Vector y(3, 1);
  y(0, 0) = 1;
  y(1, 0) = 2;
  y(2, 0) = 3;

  LinearRegression lr;
  if (lr.FitLinearModel(X, y)) {
    printf("Coefficients (B):\n");
    for (int i = 0; i < lr.B.data.size(); ++i)
      printf("  B[%d] = %.4f\n", i, lr.B.data[i]);

    Vector pvals = lr.GetAsyPvalue();
    printf("P-values:\n");
    for (int i = 0; i < pvals.Length(); ++i)
      printf("  p[%d] = %.4f\n", i, pvals[i]);

  } else {
    printf("Linear model fit failed.\n");
  }

  return 0;
}
