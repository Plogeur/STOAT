// This file is part of STOAT 0.0.1, copyright (C) 2024-2025 Matis Alias-Bagarre, Jean Monlong.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
std::vector<std::string> LMM_binary(const std::vector<std::vector<size_t>>& df,
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
std::string chi2Test(const std::vector<std::vector<size_t>>& observed) {
    
    size_t cols = observed[0].size();
    std::vector<double> col_sums(cols, 0.0);
    double row_sum0 = 0.0, row_sum1 = 0.0, total = 0.0;

    // Precompute row sums and column sums
    for (size_t j = 0; j < cols; ++j) {
        double a = observed[0][j];
        double b = observed[1][j];
        row_sum0 += a;
        row_sum1 += b;
        col_sums[j] = a + b;
        total += a + b;
    }

    if (total == 0.0) return "0.0";

    // Compute chi-squared statistic
    double chi2 = 0.0;
    for (size_t j = 0; j < cols; ++j) {
        double expected0 = (row_sum0 * col_sums[j]) / total;
        double expected1 = (row_sum1 * col_sums[j]) / total;

        double diff0 = observed[0][j] - expected0;
        double diff1 = observed[1][j] - expected1;

        if (expected0 > 0) chi2 += (diff0 * diff0) / expected0;
        if (expected1 > 0) chi2 += (diff1 * diff1) / expected1;
    }

    size_t df = (observed.size() - 1) * (observed[0].size() - 1);
    boost::math::chi_squared dist(df);
    return set_precision(boost::math::cdf(boost::math::complement(dist, chi2)));
}

// ------------------------ Fisher exact test ------------------------

std::string fastFishersExactTest(const std::vector<std::vector<size_t>>& table) {
// plink 1.9 fisher22 implementation

    // Ensure the table is 2x2
    if (table.size() != 2 || table[0].size() != 2 || table[1].size() != 2) {
        return "NA";
    }

    // Extract values from the table
    size_t m11 = table[0][0];
    size_t m12 = table[0][1];
    size_t m21 = table[1][0];
    size_t m22 = table[1][1];

    double tprob = (1 - kExactTestEpsilon2) * kExactTestBias;
    double cur_prob = tprob;
    double cprob = 0;
    size_t uii;
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
    if ((S_CAST(size_t, m11) * m22) > (S_CAST(size_t, m12) * m21)) {
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

void binary_stat_test(const std::vector<std::vector<size_t>>& df, 
    string& fastfisher_p_value, string& chi2_p_value, string& group_paths,
    string& allele_number_str, string& min_row_index_str, string& numb_colum_str, 
    string& inter_group_str, string& average_str) {

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
        size_t col_min = INT_MAX;
        for (const auto& row : df) {
            col_min = std::min(col_min, row[col]);
        }
        inter_group += col_min;
    }
    
    int average = static_cast<double>(allele_number) / numb_colum; // get 200 instead of 200.00000

    // Compute  Fisher's exact & Chi-squared test p-value
    chi2_p_value = chi2Test(df);
    fastfisher_p_value = fastFishersExactTest(df);
    group_paths = format_group_paths(df);
    allele_number_str = std::to_string(allele_number);
    min_row_index_str = std::to_string(min_row_index);
    numb_colum_str = std::to_string(numb_colum);
    inter_group_str = std::to_string(inter_group);
    average_str = std::to_string(average);
}

std::string format_group_paths(const std::vector<std::vector<size_t>>& matrix) {

    std::string result;
    size_t rows = matrix.size();

    for (size_t row = 0; row < rows; ++row) {
        result += std::to_string(matrix[row][0]) + ":" + std::to_string(matrix[row][1]);
        if (row < rows - 1) {
            result += ","; // Separate row pairs with ','
        }
    }

    return result;
}

bool create_binary_table_with_maf_check(
    std::vector<std::vector<size_t>>& df,
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, 
    const Matrix& matrix,
    const double& maf) {

    size_t length_column_headers = list_path_snarl.size();
    std::vector<size_t> g0(length_column_headers, 0);
    std::vector<size_t> g1(length_column_headers, 0);
    size_t totalSum = 0;

    for (size_t idx_g = 0; idx_g < list_path_snarl.size(); ++idx_g) {
        const std::string& path_snarl = list_path_snarl[idx_g];
        size_t number_sample = list_samples.size();

        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, number_sample * 2);

        for (int idx : idx_srr_save) {
            std::string srr = list_samples[idx / 2];

            auto it = groups.find(srr);
            if (it->second) {
                g1[idx_g] += 1;
            } else {
                g0[idx_g] += 1;
            }
            totalSum += 1;
        }
    }

    // Check MAF threshold
    for (size_t i = 0; i < length_column_headers; ++i) {
        int columnSum = g0[i] + g1[i];
        if (static_cast<double>(columnSum) / totalSum >= maf) {
            return true; // MAF threshold met
        }
    }

    df = {g0, g1};
    return false; // No column met MAF threshold
}
