#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"

using namespace std;

// Linear regression function OLS
void linear_regression(
    const std::unordered_map<std::string, std::vector<size_t>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str) {

    size_t num_samples = df.size();
    size_t max_paths = df.begin()->second.size();
    
    Eigen::MatrixXd X(num_samples, max_paths);
    X.setZero(); // Initialize matrix with zeros
    Eigen::VectorXd y(num_samples);
    
    int row = 0;
    for (const auto& [sample, paths] : df) {
        y(row) = quantitative_phenotype.at(sample);
        for (size_t col = 0; col < paths.size(); ++col) {
            X(row, col) = paths[col];
        }
        ++row;
    }
    
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;
    
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = max_paths - 1;
    int df_res = num_samples - max_paths;
    double mse = rss / df_res;  // Mean Squared Error (MSE)
    
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // Compute F-statistic
    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    boost::math::fisher_f dist(df_reg, df_res);
    double p_value = boost::math::cdf(boost::math::complement(dist, std::abs(f_stat)));

    // set precision : 4 digit
    r2_str = set_precision(r2);
    beta_str = set_precision(beta.mean());
    se_str = set_precision(se.mean());
    p_value_str = set_precision(p_value);
}

void glm_quantitative(
    const std::unordered_map<std::string, std::vector<size_t>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covar,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str) {

    size_t num_samples = df.size();
    size_t num_features = df.begin()->second.size();

    // Assuming all covariates are of the same size
    size_t num_covariates = covar.begin()->second.size();

    Eigen::MatrixXd X(num_samples, num_features + num_covariates);
    X.setZero();
    Eigen::VectorXd y(num_samples);

    int row = 0;
    for (const auto& [sample, features] : df) {
        y(row) = quantitative_phenotype.at(sample);

        // Add features (e.g., genotype/path values)
        for (size_t col = 0; col < features.size(); ++col) {
            X(row, col) = features[col];
        }

        // Add covariates
        if (covar.find(sample) != covar.end()) {
            const auto& covariate_values = covar.at(sample);
            for (size_t i = 0; i < covariate_values.size(); ++i) {
                X(row, num_features + i) = covariate_values[i];
            }
        }

        ++row;
    }

    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = num_features + num_covariates - 1;
    int df_res = num_samples - (num_features + num_covariates);
    double mse = rss / df_res;

    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    boost::math::fisher_f dist(df_reg, df_res);
    double p_value = boost::math::cdf(boost::math::complement(dist, std::abs(f_stat)));

    // set precision : 4 digits
    r2_str = set_precision(r2);
    beta_str = set_precision(beta.mean());
    se_str = set_precision(se.mean());
    p_value_str = set_precision(p_value);

}

// Function to create the quantitative table
std::pair<std::unordered_map<std::string, std::vector<size_t>>, size_t> create_quantitative_table(
    const std::vector<std::string>& list_samples, 
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    // Retrieve row headers dictionary
    size_t allele_number = 0;
    size_t length_sample = list_samples.size();
    size_t length_column = column_headers.size();

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<size_t>> genotypes(length_sample, std::vector<size_t>(length_column, 0));

    // Genotype paths
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, length_sample*2);

        for (auto idx : idx_srr_save) {
            size_t srr_idx = idx / 2;  // Adjust index to correspond to the sample index
            genotypes[srr_idx][col_idx] += 1;
            allele_number++; // inversed matrice n*m by m*n
        }
    }

    std::unordered_map<std::string, std::vector<size_t>> df;
    for (size_t i = 0; i < list_samples.size(); ++i) {
        df[list_samples[i]] = genotypes[i];
    }
    
    return {df, allele_number};
}
