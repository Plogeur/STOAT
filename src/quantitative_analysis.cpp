#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

// Student's t-distribution approximation for p-value calculation
double calculate_t_statistic(double beta, double se) {
    return beta / se;
}

// Incomplete beta function approximation for p-value calculation
double calculate_p_value(double t_statistic, int degrees_of_freedom) {
    // Simple approximation using a standard approach
    // This is a simplified p-value calculation
    double t = std::abs(t_statistic);
    double df = static_cast<double>(degrees_of_freedom);
    
    // Approximation formula for p-value
    double p = std::pow(1.0 + (t * t / df), -df/2.0);
    return 2.0 * (1.0 - p);
}

// Linear regression function
std::tuple<double, double, double> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {
    
    std::vector<double> X, Y;
    
    // Populate X and Y vectors
    for (const auto& entry : df) {
        const std::string& key = entry.first;
        auto it = quantitative_phenotype.find(key);
        if (it != quantitative_phenotype.end()) {
            for (int val : entry.second) {
                X.push_back(static_cast<double>(val));
                Y.push_back(it->second);
            }
        }
    }
    
    // Validate input data
    if (X.size() < 2 || Y.size() < 2 || X.size() != Y.size() || 
        std::adjacent_find(X.begin(), X.end(), std::not_equal_to<>()) == X.end()) {
        return std::make_tuple(0.0, 0.0, 1.0);
    }
    
    // Calculate means
    double x_mean = std::accumulate(X.begin(), X.end(), 0.0) / X.size();
    double y_mean = std::accumulate(Y.begin(), Y.end(), 0.0) / Y.size();
    
    // Calculate sum of squares
    double ssxx = 0.0, ssxy = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        double x_diff = X[i] - x_mean;
        double y_diff = Y[i] - y_mean;
        ssxx += x_diff * x_diff;
        ssxy += x_diff * y_diff;
    }
    
    // Calculate beta (slope)
    double beta = ssxy / ssxx;
    
    // Calculate standard error
    double residual_sum_sq = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        double predicted_y = x_mean + beta * (X[i] - x_mean);
        double residual = Y[i] - predicted_y;
        residual_sum_sq += residual * residual;
    }
    
    // Degrees of freedom
    int n = X.size();
    int df_ = n - 2;
    
    // Standard error calculation
    double se = std::sqrt(residual_sum_sq / df_) / std::sqrt(ssxx);
    
    // Calculate t-statistic and p-value
    double t_statistic = calculate_t_statistic(beta, se);
    double p_value = calculate_p_value(t_statistic, df_);
    
    return std::make_tuple(beta, se, p_value);
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
