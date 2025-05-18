#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cmath>

void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype) {

    size_t num_samples = df.size();
    size_t num_features = df[0].size();

    // Create matrix X with intercept
    Eigen::MatrixXd X(num_samples, num_features + 1);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept
    for (size_t i = 0; i < num_samples; ++i) {
        for (size_t j = 0; j < num_features; ++j) {
            X(i, j + 1) = df[i][j];
        }
    }

    // Response vector y
    Eigen::VectorXd y(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        y(i) = quantitative_phenotype[i];
    }

    // X^T X and X^T y
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::VectorXd Xty = X.transpose() * y;

    // Solve (X^T X) beta = X^T y using LDLT
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
    Eigen::VectorXd beta = ldlt.solve(Xty);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R-squared
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).square().sum();
    double r2 = 1.0 - (rss / tss);
    std::cout << "R²: " << r2 << std::endl;

    int df_model = X.cols() - 1;              // Exclude intercept
    int df_resid = static_cast<int>(num_samples) - static_cast<int>(X.cols());  // n - (k + 1)
    double mse = rss / df_resid;

    // Covariance matrix and standard errors
    Eigen::MatrixXd cov_matrix = ldlt.solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt();

    // Z-statistics and p-values
    std::vector<double> p_values(beta.size());
    for (int i = 0; i < beta.size(); ++i) {

        double Z = beta(i) / se(i);
        double chi2 = Z * Z;

        if (std::isnan(chi2) || chi2 < 0.0) {
            p_values[i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        boost::math::chi_squared dist(1);
        p_values[i] = boost::math::cdf(boost::math::complement(dist, chi2));

        std::cout << "beta[" << i << "]: " << beta(i)
                  << ", se: " << se(i)
                  << ", chi²: " << chi2
                  << ", p-value: " << p_values[i] << std::endl;
    }
}

int main() {
    std::vector<std::vector<double>> X = {
        {0}, {0}, {0}, {0}, {0}, {0},
        {1}, {1},
        {0.5}, {0.5}, {0.5}, {0.5},
        {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}
    };

    std::vector<double> y = {0,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0};

    linear_regression(X, y);
    return 0;
}

// g++ -std=c++17 -I/usr/local/include/eigen3 -lboost_math_c99 -o linear_regression linear_regression.cpp
