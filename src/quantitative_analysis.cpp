#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

using namespace std;

// Linear regression function
std::tuple<double, double, double> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {
    
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
    
    // Compute standard errors
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (rss / (num_samples - max_paths)) * cov_matrix.diagonal().array().sqrt().matrix();
    
    // Compute F-statistic
    int df_reg = max_paths - 1;
    int df_res = num_samples - max_paths;
    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    
    double p_value = 1.0;
    if (f_stat > 0) {
        boost::math::fisher_f dist(df_reg, df_res);
        p_value = boost::math::cdf(boost::math::complement(dist, f_stat));
    }
    
    return {beta.mean(), se.mean(), p_value};
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
