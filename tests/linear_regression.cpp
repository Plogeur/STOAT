#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
#include <boost/math/distributions/chi_squared.hpp>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <algorithm>

using namespace std;

// #define DECLARE_EIGEN_MATRIX(matRef, varName)                              \
//   Eigen::Map<Eigen::MatrixXd> varName((matRef).data.data(), (matRef).rows, \
//                                       (matRef).cols);
// #define DECLARE_EIGEN_CONST_MATRIX(matRef, varName)               \
//   Eigen::Map<const Eigen::MatrixXd> varName((matRef).data.data(), \
//                                             (matRef).rows, (matRef).cols);

// #define DECLARE_EIGEN_VECTOR(matRef, varName)               \
//   Eigen::Map<Eigen::MatrixXd> varName((matRef).data.data(), \
//                                       (matRef).data.size(), 1);

// #define DECLARE_EIGEN_CONST_VECTOR(matRef, varName)               \
//   Eigen::Map<const Eigen::MatrixXd> varName((matRef).data.data(), \
//                                             (matRef).data.size(), 1);

// // column-major
// class Matrix {
//     int rows;
//     int cols;
//     std::vector<double> data;
//     std::vector<std::string> colLabel;
// }

// class LinearRegression {
//     Vector B;     // coefficient vector
//     Matrix covB;  // coefficient covariance matrix
//     Vector pValue;     //
//     Vector residuals;  // Y - X' \hat(beta)
//     Vector predict;    // Y - X'
//                         // \hat(beta)rvtest.1110.tgzrvtest.1110.tgzrvtest.1110.tgz
//     Matrix XtXinv;     // (X'X)^ {-1}
//     double sigma2;  // \hat{\sigma^2} MLE
// }

// // use Wald statistics
// class LinearRegression {
//  public:
//   LinearRegression() : sigma2(0.){};
//   ~LinearRegression(){};

//   bool FitLinearModel(const Matrix& X,
//                       const Vector& y);  // return false if not converging

//   // alias function
//   bool Fit(Matrix& X, Matrix& y) { return this->FitLinearModel(X, y); }
//   bool Fit(const Matrix& X, const Vector& y) {
//     return this->FitLinearModel(X, y);
//   }

//   bool calculateResidualMatrix(Matrix& X, Matrix* out);
//   bool calculateHatMatrix(Matrix& X, Matrix* out);

//   Vector& GetAsyPvalue();
//   Vector& GetCovEst() { return this->B; };  // (X'X)^{-1} X'Y
//   Matrix& GetCovB() { return this->covB; };
//   Vector& GetPredicted() { return this->predict; };
//   Vector& GetResiduals() { return this->residuals; };
//   double GetSigma2() const { return this->sigma2; };
// };

// bool LinearRegression::FitLinearModel(const Matrix& X, const Vector& y) {

//   XtXinv.Dimension(X.cols, X.cols);
//   DECLARE_EIGEN_CONST_MATRIX(X, X_e);
//   DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
//   XtXinv_e = (X_e.transpose() * X_e)
//                  .llt()
//                  .solve(Eigen::MatrixXd::Identity(X_e.cols(), X_e.cols()));

//   B.Dimension(X.cols, 1);
//   DECLARE_EIGEN_VECTOR(B, B_e);
//   DECLARE_EIGEN_CONST_VECTOR(y, y_e);
//   B_e = XtXinv_e * X_e.transpose() * y_e;

//   this->predict.Dimension(X.rows, 1);
//   this->residuals.Dimension(X.rows, 1);
//   DECLARE_EIGEN_VECTOR(this->predict, predict_e);
//   DECLARE_EIGEN_VECTOR(this->residuals, resid_e);
//   predict_e = X_e * B_e;
//   resid_e = y_e - predict_e;

//   this->sigma2 = resid_e.squaredNorm() / y_e.size();
//   this->covB.Dimension(X.cols, X.cols);
//   DECLARE_EIGEN_MATRIX(this->covB, covB_e);
//   covB_e = XtXinv_e * sigma2;

//   return true;
// };

// Vector& LinearRegression::GetAsyPvalue() {
//   int numCov = B.Length();
//   pValue.Dimension(B.Length());
//   for (int i = 0; i < numCov; i++) {
//     double Zstat = B[i] / sqrt(covB(i, i));
//     Zstat *= Zstat;
//     pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
//   }
//   return (pValue);
// }

// bool LinearRegression::calculateResidualMatrix(Matrix& X, Matrix* out) {
//   if (!calculateHatMatrix(X, out)) return false;
//   Matrix& m = *out;
//   for (int i = 0; i < m.rows; ++i) {
//     for (int j = 0; j < m.cols; ++j) {
//       if (i == j) {
//         m(i, j) = 1.0 - m(i, j);
//       } else {
//         m(i, j) = -m(i, j);
//       }
//     }
//   }
//   return true;
// }

// bool LinearRegression::calculateHatMatrix(Matrix& X, Matrix* out) {

//   DECLARE_EIGEN_CONST_MATRIX(X, X_e);
//   this->XtXinv.Dimension(X.cols, X.cols);
//   DECLARE_EIGEN_MATRIX(this->XtXinv, XtXinv_e);
//   XtXinv_e = (X_e.transpose() * X_e)
//                  .llt()
//                  .solve(Eigen::MatrixXd::Identity(X_e.cols(), X_e.cols()));

//   DECLARE_EIGEN_MATRIX((*out), out_e);
//   out_e = X_e * XtXinv_e * X_e.transpose();
//   return true;
// }

void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype) {

    size_t num_samples = df.size();
    size_t max_paths = df[0].size();
    Eigen::VectorXd y(num_samples);

    for (size_t i = 0; i < num_samples; ++i) {
        y(i) = quantitative_phenotype[i];
    }

    // ------------------ WITH INTERCEPT ------------------
    // Create matrix X with intercept
    Eigen::MatrixXd X(num_samples, max_paths + 1);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept

    for (size_t i = 0; i < num_samples; ++i) {
        for (size_t j = 0; j < max_paths; ++j) {
            X(i, j + 1) = df[i][j];
        }
    }
    int df_reg = max_paths; // One additional predictor due to intercept
    int df_res = num_samples - (max_paths + 1); // Subtract all predictors + intercept
    cout << "df_res : " << df_res << endl;
    // -----------------------------------------------------

    // ----------------- WITHOUT INTERCEPT -----------------
    // Eigen::MatrixXd X(num_samples, max_paths);
    // X.setZero(); // Initialize matrix with zeros    
    // for (size_t row=0; row < num_samples; ++row) {
    //     y(row) = quantitative_phenotype[row];
    //     for (size_t col = 0; col < max_paths; ++col) {
    //         X(row, col) = df[row][col];
    //     }
    // }
    // int df_reg = max_paths - 1; // Degree of Freedom
    // int df_res = num_samples - max_paths;
    // -----------------------------------------------------

    // Coefficients beta
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    for (auto b : beta) {
        cout << "beta : " << b << endl;
    }

    // R² 
    double rss = residuals.squaredNorm();
    cout << "rss : " << rss << endl;

    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    cout << "tss : " << tss << endl;

    double r2 = 1 - (rss / tss);
    cout << "r2 : " << r2 << endl;

    double mse = rss / df_res;  // Mean Squared Error (MSE)
    cout << "mse : " << mse << endl;

    // Standard errors
    // Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse(); // WITHOUT INTER
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).ldlt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols())); // WITH INTER
    cout << "Cov diagonal: " << cov_matrix.diagonal().transpose() << endl;

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(X);
    Eigen::VectorXd beta2 = qr.solve(y);
    Eigen::MatrixXd R2 = qr.matrixR().topLeftCorner(X.cols(), X.cols());
    Eigen::MatrixXd cov_matrix2 = (R2.transpose() * R2).inverse();
    
    cout << "Cov diagonal: " << cov_matrix.diagonal().transpose() << endl;
    cout << "cov_matrix2 : " << cov_matrix2 << endl;

    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();
    Eigen::VectorXd se2 = (cov_matrix2.diagonal() * mse).array().sqrt().matrix();

    for (auto s : se) {
        cout << "se : " << s << endl;
    }

    for (auto s : se2) {
        cout << "se2 : " << s << endl;
    }
    
    // Compute F-statistic
    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    boost::math::fisher_f dist(df_reg, df_res);
    double p_value = boost::math::cdf(boost::math::complement(dist, std::abs(f_stat)));

    // t-statistics
    Eigen::VectorXd t_stats = beta.array() / se.array();
    boost::math::students_t t_dist(df_res);
    std::vector<double> p_values(beta.size());
    for (int i = 0; i < beta.size(); ++i) {
        p_values[i] = 2 * boost::math::cdf(boost::math::complement(t_dist, std::abs(t_stats[i]))); // two-tailed
        std::cout << "p-value: " << p_values[i] << std::endl;
    }
}

// python with intercept:
//                       coef    std err          t      P>|t|      [0.025      0.975]
// const               0.0692      0.052      1.327      0.186      -0.034       0.172
// >4220>4221>4223    -0.0354      0.097     -0.364      0.716      -0.227       0.156
// >4220>4222>4223     0.1047      0.123      0.854      0.394      -0.137       0.346

// python without intercept:
//                       coef    std err          t      P>|t|      [0.025      0.975]
// >4220>4221>4223     0.0338      0.097      0.349      0.727      -0.157       0.225
// >4220>4222>4223     0.1739      0.161      1.078      0.282      -0.144       0.492

// c++ with intercept :
// beta : -0.294861
// beta : 0.328683
// beta : 0.46875
// se : -nan
// se : -nan
// se : -nan

// c++ without intercept :
// beta : 0.0338218
// beta : 0.173889
// se : 0.0967754
// se : 0.161292
// p-value: 0.727095
// p-value: 0.282302

// Function to parse the feature file
void parse_feature_file(
    const std::string& feature_filename,
    std::vector<std::string>& sample_ids,
    std::vector<std::vector<double>>& X) {

    std::ifstream infile(feature_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open feature file");
    }

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token1;
        std::string token2;
        std::string sample_id;

        std::getline(ss, sample_id, '\t');
        sample_ids.push_back(sample_id);
        std::getline(ss, token1, '\t');
        std::getline(ss, token2, '\t');
        X.push_back({std::stod(token1), std::stod(token2)});
    }
}

// Function to parse the phenotype file
void parse_phenotype_file(
    const std::string& phenotype_filename,
    const std::vector<std::string>& sample_ids,
    std::vector<double>& y) {

    std::ifstream infile(phenotype_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open phenotype file");
    }

    std::unordered_map<std::string, double> phenotype_map;

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string fid, iid, pheno_str;
        std::getline(ss, fid, '\t');
        std::getline(ss, iid, '\t');
        std::getline(ss, pheno_str, '\t');
        phenotype_map[iid] = std::stod(pheno_str);
    }

    for (const auto& sample : sample_ids) {
        if (phenotype_map.find(sample) != phenotype_map.end()) {
            y.push_back(phenotype_map[sample]);
        } else {
            throw std::runtime_error("Sample ID not found in phenotype file: " + sample);
        }
    }
}

// Example usage
int main() {
    std::string feature_file = "../output/regression/4220_4223.tsv";
    std::string phenotype_file = "../data/quantitative/phenotype.tsv";

    std::vector<std::string> sample_ids;
    std::vector<std::vector<double>> X;
    std::vector<double> y;

    try {
        parse_feature_file(feature_file, sample_ids, X);
        parse_phenotype_file(phenotype_file, sample_ids, y);

        std::cout << "Parsed " << X.size() << " samples with "
                  << X[0].size() << " features.\n";
        std::cout << "Parsed " << y.size() << " y values.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    std::cout << "C++ - X:" << std::endl;
    for (size_t i = 0; i < X.size(); ++i) {
        std::cout << "Row " << i << ": ";
        for (size_t j = 0; j < X[i].size(); ++j) {
            std::cout << X[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "C++ - quantitative_phenotype[0:5]:" << std::endl;
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << y[i] << " ";
    }
    std::cout << std::endl;

    linear_regression(X, y);
    return 0;
}

// g++ -std=c++17 -I/usr/local/include/eigen3 -lboost_math_c99 -o linear_regression linear_regression.cpp
// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o linear_regression linear_regression.cpp

// Python - X:
//      const  >4220>4221>4223  >4220>4222>4223
// 0      1.0              0.0              1.0
// 1      1.0              0.0              1.0
// 2      1.0              0.0              1.0
// 3      1.0              0.0              1.0
// 4      1.0              0.0              1.0
// ..     ...              ...              ...
// 195    1.0              0.5              0.5
// 196    1.0              1.0              0.0
// 197    1.0              1.0              0.0
// 198    1.0              0.5              0.5
// 199    1.0              0.5              0.5

// [200 rows x 3 columns]
// Python - y:
// 0     -0.406483
// 1      1.008057
// 2     -1.905898
// 3     -1.591528
// 4      0.183614
//          ...   
// 195   -0.316744
// 196   -1.460248
// 197   -0.212440
// 198    1.146369
// 199    0.815490