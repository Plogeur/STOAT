#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
#include <boost/math/distributions/chi_squared.hpp>
#include <string>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <algorithm>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Linear regression function OLS with intercept
void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype) {

    size_t num_samples = df.size();
    size_t num_features = df[0].size();

    // Add intercept: X with one additional column for intercept
    Eigen::MatrixXd X(num_samples, num_features); // remove 1 column num_features + 1 (intercept) -1 (remove 1 column)
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    for (size_t row = 0; row < num_samples; ++row) {
        for (size_t col = 0; col < num_features-1; ++col) { // remove 1 column 
            X(row, col + 1) = df[row][col];
        }
    }

    // Response vector
    Eigen::VectorXd y(num_samples);
    for (size_t row = 0; row < num_samples; ++row) {
        y(row) = quantitative_phenotype[row];
    }

    // Coefficients beta
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R² 
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = X.cols() - 1;              // exclude intercept from model df
    int df_res = num_samples - X.cols();    // residual degrees of freedom
    double mse = rss / df_res;

    // Standard errors
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

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
    }
}

// Function to parse the feature file
void parse_feature_file(
    const std::string& feature_filename,
    std::vector<std::string>& sample_ids,
    std::vector<std::vector<double>>& features) {

    std::ifstream infile(feature_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open feature file");
    }

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;

        std::string sample_id;
        std::getline(ss, sample_id, '\t');
        sample_ids.push_back(sample_id);

        std::vector<double> feature_row;
        while (std::getline(ss, token, '\t')) {
            feature_row.push_back(std::stod(token));
        }

        features.push_back(feature_row);
    }
}

// Function to parse the phenotype file
void parse_phenotype_file(
    const std::string& phenotype_filename,
    const std::vector<std::string>& sample_ids,
    std::vector<double>& phenotype) {

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
            phenotype.push_back(phenotype_map[sample]);
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
    std::vector<std::vector<double>> features;
    std::vector<double> phenotype;

    try {
        parse_feature_file(feature_file, sample_ids, features);
        parse_phenotype_file(phenotype_file, sample_ids, phenotype);

        std::cout << "Parsed " << features.size() << " samples with "
                  << features[0].size() << " features.\n";
        std::cout << "Parsed " << phenotype.size() << " phenotype values.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    linear_regression(features, phenotype);
    return 0;
}

// python without intercept:
//                       coef    std err          t      P>|t|      [0.025      0.975]
// >4220>4221>4223     0.0338      0.097      0.349      0.727      -0.157       0.225
// >4220>4222>4223     0.1739      0.161      1.078      0.282      -0.144       0.492

// python with intercept + remove 1 column :
//                       coef    std err          t      P>|t|      [0.025      0.975]
// const               0.0692      0.052      1.327      0.186      -0.034       0.172
// >4220>4221>4223    -0.0354      0.097     -0.364      0.716      -0.227       0.156
// >4220>4222>4223     0.1047      0.123      0.854      0.394      -0.137       0.346

// g++ -std=c++17 -I/usr/local/include/eigen3 -lboost_math_c99 -lgsl -lgslcblas -o linear_regression linear_regression.cpp
