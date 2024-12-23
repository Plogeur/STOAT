#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

// Helper function to calculate the mean of a vector
double mean(const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// Helper function to calculate the variance of a vector
double variance(const std::vector<double>& v) {
    double m = mean(v);
    double var = 0.0;
    for (double val : v) {
        var += (val - m) * (val - m);
    }
    return var / (v.size() - 1);
}

// Helper function to calculate the covariance between two vectors
double covariance(const std::vector<double>& x, const std::vector<double>& y) {
    double mean_x = mean(x);
    double mean_y = mean(y);
    double cov = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        cov += (x[i] - mean_x) * (y[i] - mean_y);
    }
    return cov / (x.size() - 1);
}

// Linear regression function returning a tuple of p_value, standard error (se), and beta
std::tuple<double, double, double> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, float>& quantitative_phenotype) {

    // Ensure that we have matching keys in both maps
    std::vector<double> X;
    std::vector<double> Y;

    for (const auto& entry : df) {
        const std::string& key = entry.first;
        const std::vector<int>& x_values = entry.second;

        // Check if the key exists in the phenotype map
        auto it = quantitative_phenotype.find(key);
        if (it != quantitative_phenotype.end()) {
            // Convert integers to doubles and store in X
            for (int val : x_values) {
                X.push_back(static_cast<double>(val));
            }
            // Store corresponding phenotype (Y values)
            Y.push_back(static_cast<double>(it->second));
        }
    }
    
    // Print X
    std::cout << "X: ";
    for (const double& value : X) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Print Y
    std::cout << "Y: ";
    for (const double& value : Y) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    if (X.size() < 2 || Y.size() < 2 || X.size() != Y.size()) {
        throw std::invalid_argument("Data mismatch or insufficient data for linear regression.");
    }

    // Calculate the beta (slope) using the least squares method
    double beta = covariance(X, Y) / variance(X);

    // Calculate the intercept (alpha)
    double alpha = mean(Y) - beta * mean(X);

    // Calculate the sum of squared residuals (error terms)
    double sum_errors_squared = 0.0;
    for (size_t i = 0; i < X.size(); ++i) {
        double predicted_y = alpha + beta * X[i];
        double error = Y[i] - predicted_y;
        sum_errors_squared += error * error;
    }

    // Calculate the standard error (SE) of beta
    double se = std::sqrt(sum_errors_squared / (X.size() - 2)) / std::sqrt(variance(X) * (X.size() - 1));

    // Calculate t-statistic for beta and the two-tailed p-value using normal approximation
    double t_stat = beta / se;
    double p_value = 2 * (1.0 - std::erf(std::abs(t_stat) / std::sqrt(2)));  // Using the error function approximation

    return std::make_tuple(se, beta, p_value);
}

// Function to create the quantitative table
std::unordered_map<std::string, std::vector<int>> create_quantitative_table(
    const std::vector<std::string>& column_headers, 
    const std::vector<std::string>& list_samples,
    Matrix& matrix) {

    // Retrieve row headers dictionary
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();
    int length_sample = list_samples.size();
    std::vector<int> srr_save(length_sample); // replace list(range(length_sample)) 

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<int>> genotypes(length_sample, std::vector<int>(column_headers.size(), 0));

    // Fill in the matrix
    for (size_t col_idx = 0; col_idx < column_headers.size(); ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, 
                                                              matrix, length_sample*2);

        for (int idx : idx_srr_save) {
            int srr_idx = idx / 2;  // Adjust index to correspond to the sample index
            genotypes[srr_idx][col_idx] += 1;
        }
    }

    // Convert genotypes matrix to a map for easier access by sample name
    std::unordered_map<std::string, std::vector<int>> df;
    for (size_t i = 0; i < list_samples.size(); ++i) {
        df[list_samples[i]] = genotypes[i];
    }
    
    return df;
}
