#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

using namespace std;

// Linear regression function OLS
std::tuple<string, string, string, string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {
    
    if (df.empty() || quantitative_phenotype.empty()) {
        return {"Error", "Empty dataframe or phenotype", "", ""};
    }
    
    int n = df.begin()->second.size(); // Number of samples
    int p = df.size(); // Number of SNPs
    
    Eigen::MatrixXd X(n, p); // No intercept column
    Eigen::VectorXd y(n);
    
    int col_idx = 0;
    for (const auto& [snp, values] : df) {
        for (int i = 0; i < n; ++i) {
            X(i, col_idx) = values[i];
        }
        col_idx++;
    }
    
    // Fill phenotype vector
    int row_idx = 0;
    for (const auto& [sample, phenotype_value] : quantitative_phenotype) {
        y(row_idx) = phenotype_value;
        row_idx++;
    }
    
    // Compute OLS: beta = (X'X)^(-1) X'Y
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    
    // Compute residuals
    Eigen::VectorXd residuals = y - X * beta;
    double ss_res = residuals.squaredNorm();
    
    // Compute total sum of squares
    double y_mean = y.mean();
    double ss_tot = (y.array() - y_mean).square().sum();
    
    // Compute R-squared
    double r_squared = 1 - (ss_res / ss_tot);
    
    // Compute F-statistic
    double ms_model = (ss_tot - ss_res) / p;
    double ms_residual = ss_res / (n - p);
    double f_stat = ms_model / ms_residual;
    
    // Compute global p-value using Fisher's F-distribution
    boost::math::fisher_f dist(p, n - p);
    double global_p_value = 1 - boost::math::cdf(dist, f_stat);
    
    // Compute standard error
    Eigen::MatrixXd XtX_inv = (X.transpose() * X).inverse();
    Eigen::VectorXd se = XtX_inv.diagonal().array().sqrt();
    
    // Set precision to 4 digits
    string r2_str = set_precision(r_squared);
    string beta_mean_str = set_precision(beta.mean());
    string se_mean_str = set_precision(se.mean());
    string p_value_str = set_precision(global_p_value);

    return {r2_str, beta_mean_str, se_mean_str, p_value_str};
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
