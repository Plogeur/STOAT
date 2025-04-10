#include "binary_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

#include <vector>
#include <unordered_map>
#include <string>
#include <Eigen/Dense>
#include <cmath>

// Fisher's Exact Test for 2x2 contingency table
#ifndef DBL_MAX
#  define DBL_MAX 1.7976931348623157e308
#endif

#ifdef __cplusplus
#  define K_CAST(type, val) (const_cast<type>(val))
#  define R_CAST(type, val) (reinterpret_cast<type>(val))
#  define S_CAST(type, val) (static_cast<type>(val))
#endif

// 2^{-40} for now, since 2^{-44} was too small on real data
static const double kExactTestEpsilon2 = 0.0000000000009094947017729282379150390625;
static const double kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;

// ------------------------ LMM BINARY ------------------------

// LMM Model fitting function
std::vector<std::string> LMM_binary(const std::vector<std::vector<uint64_t>>& df,
    const std::unordered_map<std::string, std::vector<double>>& covariate) {

    int num_samples = df[0].size();

    // Convert df to Eigen matrix (phenotype)
    Eigen::VectorXd Y(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        Y(i) = df[1][i];  // Group 1 (affected cases)
    }

    // Covariate matrix
    int num_covariates = covariate.size();
    Eigen::MatrixXd X(num_samples, num_covariates + 1); // Include intercept

    // Set intercept
    X.col(0) = Eigen::VectorXd::Ones(num_samples);

    // Fill in covariates
    int col_idx = 1;
    for (const auto& [key, values] : covariate) {
        for (int i = 0; i < num_samples; ++i) {
            X(i, col_idx) = values[i];
        }
        col_idx++;
    }

    // Fit LMM (simplified REML method)
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
    Eigen::VectorXd residuals = Y - (X * beta);
    double sigma2 = (residuals.transpose() * residuals)(0, 0) / (num_samples - num_covariates - 1);

    // Standard errors (diagonal of covariance matrix)
    Eigen::VectorXd se = (X.transpose() * X).inverse().diagonal().array().sqrt() * std::sqrt(sigma2);

    // Compute p-value using likelihood ratio test (LRT)
    Eigen::VectorXd chi2_stat = (beta.array().square() / se.array().square());
    Eigen::VectorXd p_value = (-0.5 * chi2_stat.array()).exp(); // Approximation

    // Convert to strings for output
    return {std::to_string(beta.mean()), std::to_string(se.mean()), std::to_string(p_value.mean())};
}

// ------------------------ Chi2 test ------------------------

// Check if the observed matrix is valid (no zero rows/columns)
bool check_observed(const std::vector<std::vector<uint64_t>>& observed, uint64_t rows, uint64_t cols) {
    std::vector<int> col_sums(cols, 0);

    if (observed.size() == 0) return false;
    
    for (uint64_t i = 0; i < rows; ++i) {
        int row_sum = 0;
        for (uint64_t j = 0; j < cols; ++j) {
            row_sum += observed[i][j];
            col_sums[j] += observed[i][j];
        }
        if (row_sum <= 0) return false; // Check row sum
    }

    for (uint64_t j = 0; j < cols; ++j) {
        if (col_sums[j] <= 0) return false; // Check column sums
    }

    return true;
}

// Function to calculate the Chi-square test statistic
std::string chi2Test(const std::vector<std::vector<uint64_t>>& observed) {
    uint64_t rows = observed.size();
    uint64_t cols = observed[0].size();

    // Validate the observed matrix
    if (!check_observed(observed, rows, cols)) {
        return {"NA", "NA"};
    }

    // Compute row and column sums
    std::vector<double> row_sums(rows, 0.0);
    std::vector<double> col_sums(cols, 0.0);
    double total_sum = 0.0;

    for (uint64_t i = 0; i < rows; ++i) {
        for (uint64_t j = 0; j < cols; ++j) {
            row_sums[i] += observed[i][j];
            col_sums[j] += observed[i][j];
            total_sum += observed[i][j];
        }
    }

    // Compute expected frequencies
    std::vector<std::vector<double>> expected(rows, std::vector<double>(cols, 0.0));
    for (uint64_t i = 0; i < rows; ++i) {
        for (uint64_t j = 0; j < cols; ++j) {
            expected[i][j] = (row_sums[i] * col_sums[j]) / total_sum;
        }
    }

    // Compute chi-squared statistic
    double chi_squared_stat = 0.0;
    for (uint64_t i = 0; i < rows; ++i) {
        for (uint64_t j = 0; j < cols; ++j) {
            if (expected[i][j] > 0) { // Avoid division by zero
                double diff = observed[i][j] - expected[i][j];
                chi_squared_stat += diff * diff / expected[i][j];
            }
        }
    }

    uint64_t degrees_of_freedom = (rows - 1) * (cols - 1);

    // Compute p-value using Boost's chi-squared distribution
    boost::math::chi_squared chi_squared_dist(degrees_of_freedom);
    double p_value = boost::math::cdf(boost::math::complement(chi_squared_dist, chi_squared_stat));

    return set_precision(p_value);
}

// ------------------------ Fisher exact test ------------------------

std::string fastFishersExactTest(const std::vector<std::vector<uint64_t>>& table) {

    // Ensure the table is 2x2
    if (table.size() != 2 || table[0].size() != 2 || table[1].size() != 2) {
        return "NA";
    }

    // Extract values from the table
    uint64_t m11 = table[0][0];
    uint64_t m12 = table[0][1];
    uint64_t m21 = table[1][0];
    uint64_t m22 = table[1][1];

    double tprob = (1 - kExactTestEpsilon2) * kExactTestBias;
    double cur_prob = tprob;
    double cprob = 0;
    uint64_t uii;
    double cur11, cur12, cur21, cur22;
    double preaddp;

    if (m12 > m21) {
        uii = m12;
        m12 = m21;
        m21 = uii;
    }
    if (m11 > m22) {
        uii = m11;
        m11 = m22;
        m22 = uii;
    }
    if ((S_CAST(uint64_t, m11) * m22) > (S_CAST(uint64_t, m12) * m21)) {
        uii = m11;
        m11 = m12;
        m12 = uii;
        uii = m21;
        m21 = m22;
        m22 = uii;
    }

    cur11 = m11;
    cur12 = m12;
    cur21 = m21;
    cur22 = m22;

    while (cur12 > 0.5) {
        cur11 += 1;
        cur22 += 1;
        cur_prob *= (cur12 * cur21) / (cur11 * cur22);
        cur12 -= 1;
        cur21 -= 1;
        if (cur_prob > DBL_MAX) {
        return "0.0";
        }
        if (cur_prob < kExactTestBias) {
        tprob += cur_prob;
        break;
        }
        cprob += cur_prob;
    }

    if (cprob == 0) {
        return "1.0";
    }

    while (cur12 > 0.5) {
        cur11 += 1;
        cur22 += 1;
        cur_prob *= (cur12 * cur21) / (cur11 * cur22);
        cur12 -= 1;
        cur21 -= 1;
        preaddp = tprob;
        tprob += cur_prob;
        if (tprob <= preaddp) {
        break;
        }
    }

    if (m11) {
        cur11 = m11;
        cur12 = m12;
        cur21 = m21;
        cur22 = m22;
        cur_prob = (1 - kExactTestEpsilon2) * kExactTestBias;
        do {
        cur12 += 1;
        cur21 += 1;
        cur_prob *= (cur11 * cur22) / (cur12 * cur21);
        cur11 -= 1;
        cur22 -= 1;
        preaddp = tprob;
        tprob += cur_prob;
        if (tprob <= preaddp) {
            return set_precision(preaddp / (cprob + preaddp));
        }
        } while (cur11 > 0.5);
    }

    return set_precision(tprob / (cprob + tprob));
}

// ------------------------ Binary table & stats ------------------------

std::vector<std::string> binary_stat_test(const std::vector<std::vector<uint64_t>>& df) {

    // Compute derived statistics
    int allele_number = 0;
    int inter_group = 0;
    int numb_colum = df.empty() ? 0 : df[0].size();
    int min_row_index = INT_MAX;

    for (const auto& row : df) {
        int row_sum = std::accumulate(row.begin(), row.end(), 0);
        allele_number += row_sum;
        min_row_index = std::min(min_row_index, row_sum);
    }
    
    for (int col=0; col < numb_colum; ++col) {
        uint64_t col_min = INT_MAX;
        for (const auto& row : df) {
            col_min = std::min(col_min, row[col]);
        }
        inter_group += col_min;
    }
    
    int average = static_cast<double>(allele_number) / numb_colum; // get 200 instead of 200.00000

    // Compute  Fisher's exact & Chi-squared test p-value
    string chi2_p_value = chi2Test(df);
    string fastfisher_p_value = fastFishersExactTest(df);
    std::string group_paths = format_group_paths(df); // Placeholder for future implementation
    return {fastfisher_p_value, chi2_p_value, std::to_string(allele_number), std::to_string(min_row_index), std::to_string(numb_colum), std::to_string(inter_group), std::to_string(average), group_paths};
}

std::string format_group_paths(const std::vector<std::vector<uint64_t>>& matrix) {

    std::string result;
    uint64_t rows = matrix.size();

    for (uint64_t row = 0; row < rows; ++row) {
        result += std::to_string(matrix[row][0]) + ":" + std::to_string(matrix[row][1]);
        if (row < rows - 1) {
            result += ","; // Separate row pairs with ','
        }
    }

    return result;
}

std::vector<std::vector<uint64_t>> create_binary_table(
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, const Matrix& matrix)  {

    std::unordered_map<std::string, uint64_t> row_headers_dict = matrix.get_row_header();
    uint64_t length_column_headers = list_path_snarl.size();

    // Initialize g0 and g1 with zeros, corresponding to the length of column_headers
    std::vector<uint64_t> g0(length_column_headers, 0);
    std::vector<uint64_t> g1(length_column_headers, 0);

    // Iterate over each path_snarl in column_headers
    for (uint64_t idx_g = 0; idx_g < list_path_snarl.size(); ++idx_g) {
        const std::string& path_snarl = list_path_snarl[idx_g];
        const uint64_t number_sample = list_samples.size();
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, matrix, number_sample*2);

        // Count occurrences in g0 and g1 based on the updated idx_srr_save
        for (int idx : idx_srr_save) {
            std::string srr = list_samples[idx / 2];  // Convert index to the appropriate sample name

            // Check if sample belongs to a group and increment the respective count
            auto it = groups.find(srr);
            if (it != groups.end()) {
                if (it->second) {
                    // If true, consider as part of Group 1
                    g1[idx_g] += 1;
                } else {
                    // If false, consider as part of Group 0
                    g0[idx_g] += 1;
                }
            } else {
                throw std::runtime_error("Sample " + srr + " not found in groups.");
            }
        }
    }

    // Return the populated binary table
    return {g0, g1};
}
