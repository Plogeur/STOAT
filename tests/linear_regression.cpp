#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

// Linear regression function OLS with intercept
void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype) {

    size_t num_samples = df.size();
    size_t max_paths = df[0].size();

    Eigen::MatrixXd X(num_samples, max_paths);
    // Add intercept here 
    X.setZero(); // Initialize matrix with zeros
    Eigen::VectorXd y(num_samples);
    
    for (size_t row=0; row < num_samples; ++row) {
        y(row) = quantitative_phenotype[row];
        for (size_t col = 0; col < max_paths; ++col) {
            X(row, col) = df[row][col];
        }
    }
    
    // Coefficients beta
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    for (auto b : beta) {
        cout << "beta : " << b << endl;
    }

    // R² 
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = max_paths - 1; // Degree of Freedom
    int df_res = num_samples - max_paths;
    double mse = rss / df_res;  // Mean Squared Error (MSE)
    
    // Standard errors
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    for (auto s : se) {
        cout << "se : " << s << endl;
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
        cout << "p_values : " << p_values[i] << endl; 
    }
}

int main() {

    std::vector<std::vector<double>> X = {
    {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0}, {1, 0},
    {0, 1}, {0, 1},
    {0.5, 0.5}, {0.5, 0.5}, {0.5, 0.5}, {0.5, 0.5},
    {0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}
    };

    std::vector<double> y = {0,1,1,1,1,1,1,0,0,0,0,0,0,0,1,1,1,1,0,0};

    linear_regression(X, y);
    return 0;
}

// g++ -std=c++17 -I/usr/local/include/eigen3 -lboost_math_c99 -o linear_regression linear_regression.cpp
