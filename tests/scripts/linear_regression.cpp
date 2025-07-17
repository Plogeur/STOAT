#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <boost/math/distributions/students_t.hpp>
#include <chrono> // for benchmarking

void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_features = num_variants + 1; // +1 for intercept

    Eigen::MatrixXd X(num_samples, num_features);
    Eigen::VectorXd y(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        X(i, 0) = 1.0; // intercept
        y(i) = quantitative_phenotype[i];
        size_t col = 1;
        for (size_t j = 0; j < num_variants; ++j) {
            X(i, col++) = df[i][j];
        }
    }

    // Coefficients
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R²
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).square().sum();
    double r2 = 1 - (rss / tss);

    int df_res = static_cast<int>(num_samples - X.cols()); // residual degrees of freedom
    df_res = std::max(df_res, 1);
    double mse = rss / df_res;

    // Covariance matrix and SE
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::MatrixXd cov_matrix;
    Eigen::VectorXd se;
    cov_matrix = XtX.inverse();

    se = (cov_matrix.diagonal() * mse).array().sqrt();

    // t-stats and p-values
    Eigen::VectorXd t_stats = beta.array() / se.array();
    boost::math::students_t t_dist(df_res);

    std::vector<double> p_values;
    for (int i = 0; i < num_features; ++i) { // i = 1 avoid const p-value
        if (std::isnan(t_stats[i]) || std::isinf(t_stats[i])) {
            p_values.push_back(1.0); // Assign a high p-value for invalid t-statistics
            continue;
        }
        p_values.push_back(2 * boost::math::cdf(boost::math::complement(t_dist, std::abs(t_stats[i])))); // two-tailed
        std::cout << "p_values[" << i << "] : " << p_values[i] << std::endl;
    }

    // Print results
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Coefficients (beta):" << std::endl;
    for (int i = 0; i < num_features; ++i) {
        std::cout << "beta[" << i << "] = " << beta[i] << std::endl;
    }
    std::cout << "Standard Errors (se):" << std::endl;
    for (int i = 0; i < num_features; ++i) {
        std::cout << "se[" << i << "] = " << se[i] << std::endl;
    }
    std::cout << "R²: " << r2 << std::endl;
    std::cout << "Residual Degrees of Freedom: " << df_res << std::endl;
    std::cout << "Mean Squared Error (MSE): " << mse << std::endl;
}

// === MAIN with Example Data ===
int main() {

    std::vector<std::vector<double>> df = {
        {0.5, 0, 0.5},
        {0, 0.5, 0.5},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0.5, 0}
    };

    std::vector<double> quantitative_phenotype = {10.5, 13.0, 15.8, 19.7, 21.5};


    std::vector<std::vector<double>> covariates = {
        {1.0},
        {2.0},
        {42.0},
        {3.0},
        {2.0}
    };

    linear_regression(df, quantitative_phenotype, covariates);
    return EXIT_SUCCESS;
}

// LINUX
// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o lr linear_regression.cpp

// MACOS
// g++ -std=c++17 -I/usr/local/eigen3 -lboost_math_c99 -o lr linear_regression.cpp
