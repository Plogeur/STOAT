#include "binary_analysis.hpp"
#include "snarl_analyser.hpp"
#include "utils.hpp"

#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/math/distributions/normal.hpp>

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
static const boost::math::chi_squared chi_squared_dist(1);
boost::math::chi_squared_distribution<cpp_dec_float_50> cpp_dec_float_50_dist(1);

// ------------------------ Logistic regression ------------------------

// Standard normal cumulative distribution function
double normal_cdf(double z) {
    static const boost::math::normal_distribution<> standard_normal(0.0, 1.0);
    return boost::math::cdf(standard_normal, z);
}

// Sigmoid function
inline double sigmoid(double x) {
    return 1.0 / (1.0 + std::exp(-x));
}

// Clamp helper
inline double clamp(double x, double lo, double hi) {
    return std::max(lo, std::min(hi, x));
}

// Log-likelihood
double calculate_log_likelihood(const Eigen::VectorXd& y, const Eigen::VectorXd& p) {
    double epsilon = 1e-8;
    double ll = 0.0;
    for (int i = 0; i < y.size(); ++i) {
        double pi = clamp(p(i), epsilon, 1.0 - epsilon);
        ll += y(i) * std::log(pi) + (1 - y(i)) * std::log(1 - pi);
    }
    return ll;
}
// Logistic regression function
void logistic_regression(
    const std::vector<std::vector<double>>& variants_data,
    const std::vector<bool>& phenotype,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {

    const int max_iterations = 100;
    const double tolerance = 1e-6;
    const double l2_penalty = 1e-4;
    const double epsilon = 1e-8;

    size_t n_samples = variants_data.size();
    size_t n_variants = variants_data[0].size();

    // Add intercept column
    Eigen::MatrixXd X(n_samples, n_variants + 1);
    for (size_t i = 0; i < n_samples; ++i) {
        X(i, 0) = 1.0; // Intercept
        for (size_t j = 0; j < n_variants; ++j) {
            X(i, j + 1) = variants_data[i][j];
        }
    }

    Eigen::VectorXd y(n_samples);
    for (size_t i = 0; i < n_samples; ++i) {
        y(i) = phenotype[i] ? 1.0 : 0.0;
    }

    size_t n_params = X.cols();
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(n_params);
    Eigen::VectorXd p(n_samples);
    Eigen::VectorXd weights(n_samples);
    Eigen::VectorXd beta_old = beta;

    bool converged = false;
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd z = X * beta;
        for (int i = 0; i < n_samples; ++i) {
            p(i) = sigmoid(z(i));
            weights(i) = clamp(p(i) * (1.0 - p(i)), epsilon, 1.0);
        }

        Eigen::MatrixXd X_weighted = X;
        for (int i = 0; i < n_samples; ++i)
            X_weighted.row(i) *= std::sqrt(weights(i));

        Eigen::MatrixXd hessian = X_weighted.transpose() * X_weighted;
        hessian += l2_penalty * Eigen::MatrixXd::Identity(n_params, n_params);

        Eigen::VectorXd gradient = X.transpose() * (y - p) - l2_penalty * beta;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);
        if (ldlt.info() != Eigen::Success) return;

        Eigen::VectorXd delta = ldlt.solve(gradient);
        beta += delta;

        if ((beta - beta_old).norm() < tolerance) {
            converged = true;
            break;
        }
        beta_old = beta;
    }

    if (!converged) return;

    // Final weights
    Eigen::VectorXd z_final = X * beta;
    for (int i = 0; i < n_samples; ++i) {
        p(i) = sigmoid(z_final(i));
        weights(i) = clamp(p(i) * (1.0 - p(i)), epsilon, 1.0);
    }

    // Covariance matrix
    Eigen::MatrixXd X_weighted = X;
    for (int i = 0; i < n_samples; ++i)
        X_weighted.row(i) *= std::sqrt(weights(i));

    Eigen::MatrixXd hessian = X_weighted.transpose() * X_weighted;
    hessian += l2_penalty * Eigen::MatrixXd::Identity(n_params, n_params);
    Eigen::MatrixXd cov = hessian.inverse();
    Eigen::VectorXd se = cov.diagonal().array().sqrt();

    // p-values using Wald test
    std::vector<double> p_values(n_params);
    for (size_t i = 0; i < n_params; ++i) {
        double z_score = beta(i) / se(i);
        p_values[i] = 2.0 * (1.0 - normal_cdf(std::abs(z_score))); // Two-sided p-value
    }

    // McFadden's R²
    double ll_full = calculate_log_likelihood(y, p);
    double p_null_val = clamp(y.mean(), epsilon, 1.0 - epsilon);
    Eigen::VectorXd p_null = Eigen::VectorXd::Constant(n_samples, p_null_val);
    double ll_null = calculate_log_likelihood(y, p_null);
    double r2 = clamp(1.0 - (ll_full / ll_null), 0.0, 1.0);

    std::vector<double> p_values_adjusted = stoat_vcf::adjusted_holm(p_values);
    size_t min_index = std::distance(p_values_adjusted.begin(), std::min_element(p_values_adjusted.begin(), p_values_adjusted.end()));
    double min_p_value_adjusted = p_values_adjusted[min_index];

    // set precision : 4 digit
    r2_str = stoat_vcf::set_precision(r2);
    beta_str = stoat_vcf::set_precision(beta[min_index]);
    se_str = stoat_vcf::set_precision(se[min_index]);
    p_value_str = stoat_vcf::set_precision(min_p_value_adjusted);
}

// GLM Implementation with Iteratively Reweighted Least Squares (IRLS)
void glm_logistic_covar(
    const std::vector<std::vector<double>>& variant_data,
    const std::vector<bool>& phenotype,
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {

    const int max_iterations = 100;
    const double tolerance = 1e-6;
    const double l2_penalty = 1e-4;
    const double epsilon = 1e-8;

    size_t n_samples = variant_data.size();
    size_t n_variants = variant_data[0].size();
    size_t n_covariates = covariates[0].size();
    size_t n_features = 1 + n_variants + n_covariates; // +1 for intercept

    Eigen::MatrixXd X(n_samples, n_features);
    Eigen::VectorXd y(n_samples);

    for (size_t i = 0; i < n_samples; ++i) {
        size_t col = 0;
        X(i, col++) = 1.0; // intercept
        for (size_t j = 0; j < n_variants; ++j)
            X(i, col++) = variant_data[i][j];
        for (size_t j = 0; j < n_covariates; ++j)
            X(i, col++) = covariates[i][j];
        y(i) = phenotype[i] ? 1.0 : 0.0;
    }

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(n_features);
    Eigen::VectorXd beta_old = beta;
    Eigen::VectorXd p(n_samples);
    Eigen::VectorXd weights(n_samples);

    bool converged = false;
    for (int iter = 0; iter < max_iterations; ++iter) {
        Eigen::VectorXd z = X * beta;
        for (int i = 0; i < n_samples; ++i) {
            p(i) = sigmoid(z(i));
            weights(i) = clamp(p(i) * (1.0 - p(i)), epsilon, 1.0);
        }

        Eigen::MatrixXd X_weighted = X;
        for (int i = 0; i < n_samples; ++i)
            X_weighted.row(i) *= std::sqrt(weights(i));

        Eigen::MatrixXd hessian = X_weighted.transpose() * X_weighted;
        hessian += l2_penalty * Eigen::MatrixXd::Identity(n_features, n_features);

        Eigen::VectorXd gradient = X.transpose() * (y - p) - l2_penalty * beta;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(hessian);
        if (ldlt.info() != Eigen::Success) return;

        Eigen::VectorXd delta = ldlt.solve(gradient);
        beta += delta;

        if ((beta - beta_old).norm() < tolerance) {
            converged = true;
            break;
        }
        beta_old = beta;
    }

    if (!converged) return;

    // Final weights
    Eigen::VectorXd z_final = X * beta;
    for (int i = 0; i < n_samples; ++i) {
        p(i) = sigmoid(z_final(i));
        weights(i) = clamp(p(i) * (1.0 - p(i)), epsilon, 1.0);
    }

    // Covariance matrix
    Eigen::MatrixXd X_weighted = X;
    for (int i = 0; i < n_samples; ++i)
        X_weighted.row(i) *= std::sqrt(weights(i));

    Eigen::MatrixXd hessian = X_weighted.transpose() * X_weighted;
    hessian += l2_penalty * Eigen::MatrixXd::Identity(n_features, n_features);
    Eigen::MatrixXd cov = hessian.inverse();
    Eigen::VectorXd se = cov.diagonal().array().sqrt();

    // --- Wald Test (Normal approximation)
    std::vector<double> p_values(n_variants);
    for (size_t i = 0; i < n_variants; ++i) {
        size_t idx = 1 + i; // skip intercept
        double z_score = beta(idx) / se(idx);
        p_values[i] = 2.0 * (1.0 - normal_cdf(std::abs(z_score))); // Two-sided
    }

    // --- McFadden's R²
    double ll_full = calculate_log_likelihood(y, p);
    double p_null_val = clamp(y.mean(), epsilon, 1.0 - epsilon);
    Eigen::VectorXd p_null = Eigen::VectorXd::Constant(n_samples, p_null_val);
    double ll_null = calculate_log_likelihood(y, p_null);
    double r2 = clamp(1.0 - (ll_full / ll_null), 0.0, 1.0);

    std::vector<double> p_values_adjusted = stoat_vcf::adjusted_holm(p_values);
    size_t min_index = std::distance(p_values_adjusted.begin(), std::min_element(p_values_adjusted.begin(), p_values_adjusted.end()));
    double min_p_value_adjusted = p_values_adjusted[min_index];

    // set precision : 4 digit
    r2_str = stoat_vcf::set_precision(r2);
    beta_str = stoat_vcf::set_precision(beta[min_index]);
    se_str = stoat_vcf::set_precision(se[min_index]);
    p_value_str = stoat_vcf::set_precision(min_p_value_adjusted);
}

// ------------------------ Chi2 test ------------------------
std::string chi2_2x2(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    // Extract values from the table
    size_t a = g0[0];
    size_t b = g0[1];
    size_t c = g1[0];
    size_t d = g1[1];

    int64_t row1 = a + b;
    int64_t row2 = c + d;
    int64_t col1 = a + c;
    int64_t col2 = b + d;
    int64_t total = row1 + row2;

    if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) return "0.0";

    double expected_a = (double)(row1) * (col1) / total;
    double expected_b = (double)(row1) * (col2) / total;
    double expected_c = (double)(col1) * (row2) / total;
    double expected_d = (double)(col2) * (row2) / total;

    if (expected_a == 0 || expected_b == 0 || expected_c == 0 || expected_d == 0)
        return stoat_vcf::set_precision(std::numeric_limits<double>::max());

    double chi2_stat = 0;
    chi2_stat += std::pow((double)a - expected_a, 2) / expected_a;
    chi2_stat += std::pow((double)b - expected_b, 2) / expected_b;
    chi2_stat += std::pow((double)c - expected_c, 2) / expected_c;
    chi2_stat += std::pow((double)d - expected_d, 2) / expected_d;

    if (chi2_stat > 85.0) {
        cpp_dec_float_50 chi2_stat_float_50 = chi2_stat;
        cpp_dec_float_50 pval = 1.0 - boost::math::cdf(cpp_dec_float_50_dist, chi2_stat_float_50);
        return stoat_vcf::set_precision_float_50(pval.convert_to<double>());
    }
    return stoat_vcf::set_precision(1.0 - boost::math::cdf(chi_squared_dist, chi2_stat));
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

    size_t df = cols - 1;
    if (chi2 > 85.0) { // avoiding case 0.000+00 precision
        cpp_dec_float_50 chi2_stat_float_50 = chi2;
        boost::math::chi_squared_distribution<cpp_dec_float_50> cpp_dec_float_50_dist_2xN(df);
        cpp_dec_float_50 p_value = 1.0 - boost::math::cdf(cpp_dec_float_50_dist_2xN, chi2_stat_float_50);
        return stoat_vcf::set_precision_float_50(p_value);
    }

    boost::math::chi_squared dist_2xN(df);
    double pvalue = 1.0 - boost::math::cdf(dist_2xN, chi2);
    return stoat_vcf::set_precision(pvalue);
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

    // Ensure we are left of the distribution center, m11 <= m22, and m12 <= m21.
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
            return stoat_vcf::set_precision(preaddp / (cprob + preaddp));
        }
        } while (cur11 > 0.5);
    }

    return stoat_vcf::set_precision(tprob / (cprob + tprob));
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
    const EdgeBySampleMatrix& matrix) {

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
