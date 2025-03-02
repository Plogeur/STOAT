#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

// Function to perform linear regression using Boost
std::tuple<double, double, std::string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {
    
    std::vector<double> X, Y;
    
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
    
    if (X.size() < 2 || Y.size() < 2 || X.size() != Y.size()) {
        return std::make_tuple(0.0, 0.0, "NA");
    }
    
    auto [beta, alpha] = boost::math::statistics::simple_ordinary_least_squares(X, Y);
    
    // Compute residuals and standard error
    double sum_errors_squared = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        double predicted_y = alpha + beta * X[i];
        double error = Y[i] - predicted_y;
        sum_errors_squared += error * error;
    }
    
    double variance_x = boost::math::statistics::variance(X);
    double se = std::sqrt(sum_errors_squared / (X.size() - 2)) / std::sqrt(variance_x * (X.size() - 1));
    
    // Compute t-statistic and p-value
    double t_stat = beta / se;
    boost::math::students_t dist(X.size() - 2);
    double p_value = 2 * (1 - boost::math::cdf(dist, std::abs(t_stat)));
    
    // Format output values
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(4) << p_value;
    
    return std::make_tuple(se, beta, ss.str());
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
        // std::cout << "path_snarl : " << path_snarl << std::endl;
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, 
                                                              matrix, length_sample*2);

        for (auto idx : idx_srr_save) {
            size_t srr_idx = idx / 2;  // Adjust index to correspond to the sample index
            // std::cout << "idx : " << idx << std::endl;
            // std::cout << "srr_idx : " << srr_idx << std::endl;
            genotypes[srr_idx][col_idx] += 1;
        }
        // std::cout << std::endl;
    }

    std::unordered_map<std::string, std::vector<int>> df;
    for (size_t i = 0; i < list_samples.size(); ++i) {
        df[list_samples[i]] = genotypes[i];
    }
    
    return df;
}
