#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

using namespace std;

// Function to fit a Linear Mixed Model (LMM)
std::tuple<std::string, std::string, std::string> LMM_quantitative(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covariate) {

    size_t num_samples = df.size();
    size_t num_covariates = covariate.begin()->second.size();

    Eigen::VectorXd y(num_samples);  // Phenotype vector
    Eigen::MatrixXd X(num_samples, num_covariates + 1);  // Design matrix (including intercept)
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    
    int row = 0;
    for (const auto& [sample, paths] : df) {
        y(row) = quantitative_phenotype.at(sample);

        for (size_t col = 0; col < num_covariates; ++col) {
            X(row, col + 1) = covariate.at(sample)[col];
        }
        ++row;
    }

    // Compute the LMM parameters using OLS as an approximation
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;
    
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = num_covariates;
    int df_res = num_samples - num_covariates - 1;
    double mse = rss / df_res;

    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // Compute p-value for the model fit using Chi-Square test
    double chi2_stat = (beta.transpose() * X.transpose() * X * beta)(0, 0) / mse;
    boost::math::chi_squared dist(df_reg);
    double p_value = boost::math::cdf(boost::math::complement(dist, chi2_stat));

    // Convert to strings with precision
    std::string r2_str = set_precision(r2);
    std::string beta_mean_str = set_precision(beta.mean());
    std::string p_value_str = set_precision(p_value);

    return {r2_str, beta_mean_str, p_value_str};
}

// Linear regression function OLS
std::tuple<string, string, string, string, string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype, 
    const size_t& total_snarl) {
    
    size_t num_samples = df.size();
    size_t max_paths = 0;
    for (const auto& [_, paths] : df) {
        max_paths = std::max(max_paths, paths.size());
    }
    
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
    
    long double p_value = 1.0f;
    if (f_stat > 0) {
        boost::math::fisher_f dist(df_reg, df_res);
        p_value = boost::math::cdf(boost::math::complement(dist, f_stat));
    }

    // set precision : 4 digit
    string r2_str = set_precision(r2);
    string bete_mean_str = set_precision(beta.mean());
    string se_mean_str = set_precision(se.mean());
    string p_value_str = set_precision(p_value);
    return {r2_str, bete_mean_str, se_mean_str, p_value_str};
}

// Function to create the quantitative table
std::pair<std::unordered_map<std::string, std::vector<int>>, size_t> create_quantitative_table(
    const std::vector<std::string>& list_samples, 
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    // Retrieve row headers dictionary
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();
    size_t allele_number = 0;
    size_t length_sample = list_samples.size();
    size_t length_column = column_headers.size();
    std::vector<int> srr_save(length_sample);

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<int>> genotypes(length_sample, std::vector<int>(length_column, 0));

    // Genotype paths
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, 
                                                              matrix, length_sample*2);

        for (auto idx : idx_srr_save) {
            size_t srr_idx = idx / 2;  // Adjust index to correspond to the sample index
            genotypes[srr_idx][col_idx] += 1;
            allele_number++;
        }
    }

    std::unordered_map<std::string, std::vector<int>> df;
    for (size_t i = 0; i < list_samples.size(); ++i) {
        df[list_samples[i]] = genotypes[i];
    }
    
    return {df, allele_number};
}
