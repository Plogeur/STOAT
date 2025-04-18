#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <boost/math/distributions/fisher_f.hpp>

using namespace std;
using namespace Eigen;

void linear_regression(const MatrixXd& X, const VectorXd& y) {
    // Add intercept (bias) term
    MatrixXd X_b(X.rows(), X.cols() + 1);
    X_b << VectorXd::Ones(X.rows()), X;

    // Compute coefficients: (X'X)^(-1) X'y
    VectorXd coefficients = (X_b.transpose() * X_b).ldlt().solve(X_b.transpose() * y);

    VectorXd y_pred = X_b * coefficients;
    VectorXd residuals = y - y_pred;

    double ss_res = residuals.squaredNorm();  // RSS
    double ss_tot = (y.array() - y.mean()).square().sum();  // TSS
    double r_squared = 1 - ss_res / ss_tot;

    int n = X_b.rows();
    int k = X_b.cols() - 1;

    double msr = (ss_tot - ss_res) / k;            // Model mean square
    double mse = ss_res / (n - k - 1);             // Error mean square
    double f_statistic = msr / mse;

    // p-value from F-distribution
    boost::math::fisher_f dist(k, n - k - 1);
    double p_value = 1 - boost::math::cdf(dist, f_statistic);

    // Output
    cout << "Regression coefficients:\n" << coefficients << "\n\n";
    cout << "R-squared: " << r_squared << endl;
    cout << "F-statistic: " << f_statistic << endl;
    cout << "p-value of F-statistic: " << p_value << endl;
}

int main() {
    // Input data: samples with features
    map<string, vector<int>> data = {
        {"sample1", {1, 0, 1}},
        {"sample2", {0, 1, 1}},
        {"sample3", {1, 1, 0}},
        {"sample4", {0, 0, 1}},
        {"sample5", {1, 1, 1}}
    };

    // Phenotype target values
    map<string, double> pheno = {
        {"sample1", 2.5},
        {"sample2", 1.8},
        {"sample3", 3.0},
        {"sample4", 1.2},
        {"sample5", 3.5}
    };

    int num_samples = data.size();
    int num_features = data.begin()->second.size();

    MatrixXd X(num_samples, num_features);
    VectorXd y(num_samples);

    int row = 0;
    for (const auto& [sample, features] : data) {
        for (int col = 0; col < num_features; ++col) {
            X(row, col) = features[col];
        }
        y(row) = pheno[sample];
        row++;
    }

    linear_regression(X, y);

    return 0;
}
