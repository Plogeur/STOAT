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

using boost::multiprecision::cpp_dec_float_50;
using boost::math::chi_squared_distribution;

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

// ------------------------ Logistic regression + covariate test ------------------------

double sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

// Compute mean-centered R² (McFadden's pseudo R²)
double compute_r2(const Eigen::VectorXd& y, const Eigen::VectorXd& p_null, const Eigen::VectorXd& p_full) {
    double ll_null = (y.array() * p_null.array().log() + (1.0 - y.array()) * (1.0 - p_null.array()).log()).sum();
    double ll_full = (y.array() * p_full.array().log() + (1.0 - y.array()) * (1.0 - p_full.array()).log()).sum();
    return 1.0 - (ll_full / ll_null);
}

void logistic_regression(
// Simple Gaussian elimination to solve Ax = b
std::vector<double> solve_linear_system(std::vector<std::vector<double>> A, std::vector<double> b) {
    size_t n = A.size();

    for (size_t i = 0; i < n; ++i) {
        // Pivot
        size_t max_row = i;
        for (size_t k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                max_row = k;
            }
        }
        std::swap(A[i], A[max_row]);
        std::swap(b[i], b[max_row]);

        // Eliminate
        for (size_t k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (size_t j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

// Logistic regression
void logistic_regression_covar(
    const std::vector<std::vector<size_t>>& variant_data,
    const std::vector<bool>& phenotype,
    const std::vector<double>& covariates,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str) {

    const std::size_t n = phenotype.size();
    const std::size_t num_paths = variant_data[0].size();
    const std::size_t num_covs = covariates.size();

    // Convert phenotype to Eigen vector
    Eigen::VectorXd y(n);
    for (std::size_t i = 0; i < n; ++i) {
        y(i) = phenotype[i] ? 1.0 : 0.0;
    }

    // Sum alleles over all paths per sample
    Eigen::VectorXd snp(n);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < num_paths; ++j) {
            sum += static_cast<double>(variant_data[i][j]);
        }
        snp(i) = sum;
    }

    // Build design matrix X: intercept + covariates + SNP
    Eigen::MatrixXd X(n, num_covs + 2);
    X.col(0) = Eigen::VectorXd::Ones(n); // intercept

    std::size_t col = 1;
    for (const auto& kv : covariates) {
        const std::vector<double>& values = kv.second;
        for (std::size_t i = 0; i < n; ++i) {
            X(i, col) = values[i];
        }
        ++col;
    }

    X.col(col) = snp; // last column is SNP

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(X.cols());
    const int max_iter = 25;
    const double tol = 1e-6;

    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::VectorXd z = X * beta;
        Eigen::VectorXd p = z.unaryExpr([](double val) { return sigmoid(val); });
        Eigen::VectorXd W_diag = p.array() * (1.0 - p.array());

        // Construct diagonal weight matrix
        Eigen::MatrixXd W = W_diag.asDiagonal();

        // Compute gradient and Hessian
        Eigen::VectorXd grad = X.transpose() * (y - p);
        Eigen::MatrixXd H = X.transpose() * W * X;

        // Newton-Raphson update step
        Eigen::VectorXd delta = H.ldlt().solve(grad);
        beta += delta;

        if (delta.norm() < tol) break;
    }

    // Extract SNP coefficient
    double snp_beta = beta(beta.size() - 1);

    // Compute standard error from inverse Hessian
    Eigen::VectorXd z_final = X * beta;
    Eigen::VectorXd p_final = z_final.unaryExpr([](double val) { return sigmoid(val); });
    Eigen::VectorXd W_diag_final = p_final.array() * (1.0 - p_final.array());
    Eigen::MatrixXd W_final = W_diag_final.asDiagonal();
    Eigen::MatrixXd H_final = X.transpose() * W_final * X;

    double se = std::sqrt(H_final.inverse()(beta.size() - 1, beta.size() - 1));

    // Compute z-score and p-value
    double z_stat = snp_beta / se;
    double pval = 2.0 * (1.0 - std::erf(std::abs(z_stat) / std::sqrt(2.0)));

    // Compute pseudo-R²
    Eigen::VectorXd p_null = Eigen::VectorXd::Constant(n, sigmoid(beta(0))); // intercept-only
    Eigen::VectorXd p_full = p_final;

    double r2 = compute_r2(y, p_null, p_full);

    // Set output strings
    cout << "pval: " << pval << endl;
    p_value_str = std::to_string(pval);
    beta_str = std::to_string(snp_beta);
    se_str = std::to_string(se);
    r2_str = std::to_string(r2);
}

// ------------------------ Chi2 test ------------------------

std::string chi2_2x2(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    int a = g0[0];
    int b = g0[1];
    int c = g1[0];
    int d = g1[1];

    int row1 = a + b;
    int row2 = c + d;
    int col1 = a + c;
    int col2 = b + d;
    int total = row1 + row2;

    if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) {
        return "NA";
    }

    double numerator = static_cast<double>(a * d - b * c);
    numerator = std::abs(numerator) - 0.5 * total;
    numerator = std::max(0.0, numerator);
    numerator *= numerator;
    double denominator = static_cast<double>(row1 * row2 * col1 * col2) / total;
    long double chi2_stat = numerator / denominator;

    if (chi2_stat > 85.0) {
        cpp_dec_float_50 chi2_stat_float_50 = chi2_stat;
        chi_squared_distribution<cpp_dec_float_50> dist(1);
        cpp_dec_float_50 p_value = 1.0 - boost::math::cdf(dist, chi2_stat_float_50);
        return set_precision_chi2(p_value);
    }

    boost::math::chi_squared dist(1);
    long double p_value = 1.0 - boost::math::cdf(dist, chi2_stat);
    return set_precision(p_value);
}

// Check if the observed matrix is valid (no zero rows/columns)
std::string chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    size_t cols = g0.size();
    std::vector<size_t> col_totals(cols);
    size_t total = 0;
    size_t row_total_0 = 0;
    size_t row_total_1 = 0;

    for (size_t i = 0; i < cols; ++i) {
        col_totals[i] = g0[i] + g1[i];
        total += col_totals[i];
        row_total_0 += g0[i];
        row_total_1 += g1[i];
    }

    if (total == 0)
        return "NA";
    if (row_total_0 == 0 || row_total_1 == 0)
        return "NA";
    if (std::any_of(col_totals.begin(), col_totals.end(), [](int x){ return x == 0; }))
        return "NA";

    // Compute chi-squared
    double chi2 = 0.0;
    for (size_t i = 0; i < cols; ++i) {
        double expected_0 = static_cast<double>(row_total_0) * col_totals[i] / total;
        double expected_1 = static_cast<double>(row_total_1) * col_totals[i] / total;

        chi2 += (g0[i] - expected_0) * (g0[i] - expected_0) / expected_0;
        chi2 += (g1[i] - expected_1) * (g1[i] - expected_1) / expected_1;
    }

    if (chi2 > 85.0) { // avoiding case 0.000+00 precision
        cpp_dec_float_50 chi2_stat_float_50 = chi2;
        chi_squared_distribution<cpp_dec_float_50> dist(1);
        cpp_dec_float_50 p_value = 1.0 - boost::math::cdf(dist, chi2_stat_float_50);
        return set_precision_chi2(p_value);
    }

    size_t df = cols - 1;
    boost::math::chi_squared dist(df);
    double pvalue = 1.0 - boost::math::cdf(dist, chi2);
    return set_precision(pvalue);
}

// ------------------------ Fisher exact test ------------------------

std::string fastFishersExactTest(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {
// plink 1.9 fisher22 implementation

    // Extract values from the table
    size_t m11 = g0[0];
    size_t m12 = g0[1];
    size_t m21 = g1[0];
    size_t m22 = g1[1];

    // Check for any full-zero row or column
    if ((m11 | m12) == 0 || (m21 | m22) == 0 || (m11 | m21) == 0 || (m12 | m22) == 0) {
        return "NA";
    }
    
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
        return "1.0000";
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

void binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    string& fastfisher_p_value, string& chi2_p_value, string& group_paths,
    string& allele_number_str, string& min_row_index_str, string& numb_colum_str, 
    string& inter_group_str, string& average_str) {

    // Compute derived statistics
    int allele_number = 0;
    int inter_group = 0;
    int numb_colum = g0.size();
    int min_row_index = INT_MAX;

    for (size_t i = 0; i < g0.size(); ++i) {
        int row_sum = static_cast<int>(g0[i] + g1[i]);
        allele_number += row_sum;
        min_row_index = std::min(min_row_index, row_sum);
    }

    for (int col=0; col < numb_colum; ++col) {
        size_t col_min = INT_MAX;
        col_min = std::min(col_min, g0[col]);
        col_min = std::min(col_min, g1[col]);
        inter_group += col_min;
    }
    
    int average = static_cast<double>(allele_number) / numb_colum; // get 200 instead of 200.00000

    // Compute  Fisher's exact & Chi-squared test p-value
    if (g0.size() == 2) {
        chi2_p_value = chi2_2x2(g0, g1);
        fastfisher_p_value = fastFishersExactTest(g0, g1);
    } else {
        chi2_p_value = chi2_2xN(g0, g1);
    }
    group_paths = format_group_paths(g0, g1);
    allele_number_str = std::to_string(allele_number);
    min_row_index_str = std::to_string(min_row_index);
    numb_colum_str = std::to_string(numb_colum);
    inter_group_str = std::to_string(inter_group);
    average_str = std::to_string(average);
}

std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    std::string result;
    size_t numb_col = g0.size();
    for (size_t index_col = 0; index_col < numb_col; ++index_col) {
        result += std::to_string(g0[index_col]) + ":" + std::to_string(g1[index_col]);
        if (index_col < numb_col - 1) {
            result += ","; // Separate row pairs with ','
        }
    }
    return result;
}

size_t create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const vector<bool>& binary_phenotype, 
    const std::vector<std::string>& list_path_snarl, 
    const size_t& number_paths,
    const size_t& number_samples,
    const Matrix& matrix) {

    size_t total_sum = 0;
    for (size_t idx_g = 0; idx_g < number_paths; ++idx_g) {
        const std::string& path_snarl = list_path_snarl[idx_g];

        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        std::vector<size_t> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, number_samples * 2);

        for (size_t idx : idx_srr_save) {
            bool group = binary_phenotype[idx / 2];
            if (group) {
                g1[idx_g] += 1;
            } else {
                g0[idx_g] += 1;
            }
            total_sum++;
        }
    }
    return total_sum;
}

bool check_MAF_threshold(
    const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    const size_t& totalSum, const size_t& length_column_headers, 
    const double& maf) {

    // Check MAF threshold
    for (size_t i = 0; i < length_column_headers; ++i) {
        int columnSum = g0[i] + g1[i];
        if (static_cast<double>(columnSum) / totalSum >= maf) {
            return true; // MAF threshold met
        }
    }
    return false; // No column met MAF threshold
}