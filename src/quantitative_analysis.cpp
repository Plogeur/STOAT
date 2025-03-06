#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

// Linear regression function
std::tuple<double, double, double> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {
        
    size_t num_samples = df.size();
    size_t max_paths = 0;
    for (const auto& [sample, paths] : df) {
        if (paths.size() > max_paths) {
            max_paths = paths.size();
        }
    }
    
    MatrixXd X(num_samples, max_paths);
    X.setZero(); // Initialize matrix with zeros
    VectorXd y(num_samples);
    
    int row = 0;
    for (const auto& [sample, paths] : df) {
        y(row) = quantitative_phenotype.at(sample);
        for (size_t col = 0; col < paths.size(); ++col) {
            X(row, col) = paths[col];
        }
        row++;
    }
    
    VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    
    VectorXd y_pred = X * beta;
    VectorXd residuals = y - y_pred;
    
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);
    
    // Compute standard errors
    MatrixXd cov_matrix = (X.transpose() * X).inverse();
    VectorXd se = residuals.squaredNorm() / (num_samples - max_paths) * cov_matrix.diagonal().array().sqrt().matrix();
    
    // Compute F-statistic
    int df_reg = max_paths - 1;
    int df_res = num_samples - max_paths;
    if (df_res <= 0) {
        std::cerr << "Error: Degrees of freedom for residuals must be positive." << std::endl;
        return;
    }
    
    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    
    // Compute p-value using an approximation (F-distribution p-value calculation can be done using statistical libraries)
    double p_value = std::exp(-0.5 * f_stat);
    return std::make_tuple(beta.mean(), se.mean(), p_value);
}

// Function to create the quantitative table
std::unordered_map<std::string, std::vector<int>> create_quantitative_table(
    const std::vector<std::string>& list_samples, 
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    // Retrieve row headers dictionary
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();
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
        }
    }

    std::unordered_map<std::string, std::vector<int>> df;
    for (size_t i = 0; i < list_samples.size(); ++i) {
        df[list_samples[i]] = genotypes[i];
    }
    
    return df;
}
