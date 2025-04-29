#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"

using namespace std;

// Linear regression function OLS
void linear_regression(
    const std::vector<std::vector<size_t>>& df,
    const std::vector<double>& quantitative_phenotype,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str) {

    size_t num_samples = df.size();
    size_t max_paths = df[0].size();

    Eigen::MatrixXd X(num_samples, max_paths);
    X.setZero(); // Initialize matrix with zeros
    Eigen::VectorXd y(num_samples);
    
    for (size_t row=0; row < num_samples; ++row) {
        y(row) = quantitative_phenotype[row];
        for (size_t col = 0; col < max_paths; ++col) {
            X(row, col) = df[row][col];
        }
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
    const std::vector<std::vector<size_t>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covar,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str) 
{
    size_t num_samples = df.size();
    size_t num_features = df[0].size();  // Assumes all rows have the same number of features
    size_t num_covariates = covar.begin()->second.size();  // Assumes covariate vectors are same length

    Eigen::MatrixXd X(num_samples, num_features + num_covariates);
    Eigen::VectorXd y(num_samples);

    for (size_t row = 0; row < num_samples; ++row) {
        y(row) = quantitative_phenotype[row];

        // Add features
        for (size_t col = 0; col < num_features; ++col) {
            X(row, col) = static_cast<double>(df[row][col]);
        }

        // Add covariates (assuming covariates are indexed in the same order as samples)
        size_t i = 0;
        for (const auto& [_, covariate_vector] : covar) {
            X(row, num_features + i) = covariate_vector[row];
            ++i;
        }
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
std::pair<std::vector<std::vector<size_t>>, size_t> create_quantitative_table(
    const size_t& length_sample,
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    // Retrieve row headers dictionary
    size_t allele_number = 0;
    size_t length_column = column_headers.size();

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<size_t>> genotypes(length_sample, std::vector<size_t>(length_column, 0));

    // Genotype paths
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<size_t> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, length_sample*2);

        for (size_t idx : idx_srr_save) {
            size_t srr_idx = idx / 2;  // Adjust index to correspond to the sample index
            genotypes[srr_idx][col_idx] += 1;
            allele_number++;
        }
    }
    
    return {genotypes, allele_number};
}
