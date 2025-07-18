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
    size_t num_covariates = 0;
    size_t num_features = num_variants + 1; // +1 for intercept

    if (!covar.empty()) {
        num_covariates = covar[0].size();
        num_features = num_variants + num_covariates + 1; // +1 for intercept
    }

    Eigen::MatrixXd X(num_samples, num_features);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    Eigen::VectorXd y(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        y(i) = quantitative_phenotype[i];
        size_t col = 1;
        for (size_t j = 0; j < num_variants; ++j) {
            X(i, col++) = df[i][j];
        }
        for (size_t j = 0; j < num_covariates; ++j) {
            X(i, col++) = covar[i][j];
        }
    }
    
    // Coefficients beta
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R²
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_res = (num_samples - X.cols() + 1); // residual degrees of freedom
    df_res = std::max(df_res, 1); // Ensure df_res is at least 1 to avoid division by zero
    double mse = rss / df_res;

    // Standard errors
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();    
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // change cov_matrix calcul if X.transpose() * X might be ill-conditioned or nearly singular
    if (se.hasNaN()) {
        Eigen::MatrixXd XtX = X.transpose() * X;
        Eigen::MatrixXd cov_matrix = XtX.ldlt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
        se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();
    }

    // t-statistics
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
        {0},
        {1},
        {0}
    };

    std::vector<double> quantitative_phenotype = {2.0, 4.0, 6.0};

    std::vector<std::vector<double>> covariates = {
    };

    linear_regression(df, quantitative_phenotype, covariates);
    return EXIT_SUCCESS;
}

// LINUX
// g++ -std=c++17 -I/usr/include/eigen3 -lboost_math_c99 -o lr linear_regression.cpp

// MACOS
// g++ -std=c++17  -I/usr/local/include/eigen3 -lboost_math_c99 -o lr linear_regression.cpp

// ==============================================================================
// Dep. Variable:                      y   R-squared:                       0.000
// Model:                            OLS   Adj. R-squared:                 -1.000
// Method:                 Least Squares   F-statistic:                     0.000
// Date:                Fri, 18 Jul 2025   Prob (F-statistic):               1.00
// Time:                        11:15:37   Log-Likelihood:                -5.7281
// No. Observations:                   3   AIC:                             15.46
// Df Residuals:                       1   BIC:                             13.65
// Df Model:                           1                                         
// Covariance Type:            nonrobust                                         
// ==============================================================================
//                  coef    std err          t      P>|t|      [0.025      0.975]
// ------------------------------------------------------------------------------
// const          4.0000      2.000      2.000      0.295     -21.412      29.412
// x1          1.332e-15      3.464   3.85e-16      1.000     -44.016      44.016
// ==============================================================================
// Omnibus:                          nan   Durbin-Watson:                   1.000
// Prob(Omnibus):                    nan   Jarque-Bera (JB):                0.281
// Skew:                           0.000   Prob(JB):                        0.869
// Kurtosis:                       1.500   Cond. No.                         2.41
// ==============================================================================