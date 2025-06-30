#include "stats_test.hpp"

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

// GLM Implementation with Iteratively Reweighted Least Squares (IRLS)
void logistic_regression(
    const std::vector<std::vector<double>>& variant_data,
    const std::vector<bool>& phenotype,
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, 
    std::string& beta_str, 
    std::string& se_str, 
    std::string& r2_str) {

    const int max_iterations = 100;
    const double tolerance = 1e-6;
    const double l2_penalty = 1e-4;
    const double epsilon = 1e-8;

    size_t n_samples = variant_data.size();
    size_t n_variants = variant_data[0].size();
    size_t n_covariates = 0;
    size_t n_features = n_variants + 1; // +1 for intercept
    
    if (!covariates.empty()) {
        size_t n_covariates = covariates[0].size();
        size_t n_features =  n_variants + n_covariates + 1; // +1 for intercept
    }

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
std::string chi2_2x2(const size_t& a, const size_t& b, const size_t& c, const size_t& d) {

    int64_t row1 = a + b;
    int64_t row2 = c + d;
    int64_t col1 = a + c;
    int64_t col2 = b + d;
    int64_t total = row1 + row2;

    if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) return "NA";

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

// Fisher's exact test for a 2x2 contingency table
// m11, m12, m21, m22 are the counts in the table
// Returns the p-value as a std::string with 4 decimal places
std::string fastFishersExactTest(size_t m11, size_t m12,
                                 size_t m21, size_t m22) {
    
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

// ------------------------ Linear regression ------------------------

// Linear regression function OLS with intercept + covariate
void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar,
    std::string& p_value_str, 
    std::string& beta_str, 
    std::string& se_str, 
    std::string& r2_str) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = 0;
    size_t num_features = num_samples + 1; // +1 for intercept

    if (!covar.empty()) {
        size_t num_covariates = covar[0].size();
        size_t num_features = num_variants + num_covariates + 1; // +1 for intercept
    }

    Eigen::MatrixXd X(num_samples, num_features);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    Eigen::VectorXd y(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        y(i) = quantitative_phenotype[i];
        size_t col = 1;
        for (size_t j = 0; j < num_variants; ++j) {
            X(i, col++) = df[i][j];
        }
        for (size_t j = 0; j < num_covariates; ++j) {
            X(i, col++) = covar[i][j];
        }
    }
    
    // Coefficients beta
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R²
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_res = (num_samples - X.cols() + 1); // residual degrees of freedom
    df_res = std::max(df_res, 1); // Ensure df_res is at least 1 to avoid division by zero
    double mse = rss / df_res;

    // Standard errors
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();    
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // change cov_matrix calcul if X.transpose() * X might be ill-conditioned or nearly singular
    if (se.hasNaN()) {
        Eigen::MatrixXd XtX = X.transpose() * X;
        Eigen::MatrixXd cov_matrix = XtX.ldlt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
        se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();
        // std::cerr << "Warning: se is nan" << std::endl;
    }

    // t-statistics
    Eigen::VectorXd t_stats = beta.array() / se.array();
    boost::math::students_t t_dist(df_res);
 
    std::vector<double> p_values;
    for (int i = 1; i < num_features; ++i) { // i = 1 avoid const p-value
        if (std::isnan(t_stats[i]) || std::isinf(t_stats[i])) {
            p_values.push_back(1.0); // Assign a high p-value for invalid t-statistics
            continue;
        }
        p_values.push_back(2 * boost::math::cdf(boost::math::complement(t_dist, std::abs(t_stats[i])))); // two-tailed
    }

    double p_value_adjusted = 0;
    double beta_adjusted = 0;
    double se_adjusted = 0;

    if (p_values.size() > 1) {
        std::vector<double> p_values_adjusted = stoat_vcf::adjusted_holm(p_values);
        size_t min_index = std::distance(p_values_adjusted.begin(), std::min_element(p_values_adjusted.begin(), p_values_adjusted.end()));
        p_value_adjusted = p_values_adjusted[min_index];
        beta_adjusted = beta[min_index+1];
        se_adjusted = se[min_index+1];

    } else {
        p_value_adjusted = p_values[0];
        beta_adjusted = beta[0];
        se_adjusted = se[0];
    }

    // set precision : 4 digit
    r2_str = stoat_vcf::set_precision(r2);
    beta_str = stoat_vcf::set_precision(beta_adjusted);
    se_str = stoat_vcf::set_precision(se_adjusted);
    p_value_str = stoat_vcf::set_precision(p_value_adjusted);
}

// Eigen::MatrixXd MatrixtoEigenMatrix(const std::vector<std::vector<double>>& matrix) {
//     size_t N = matrix.size();
//     Eigen::MatrixXd M(N, N);
//     for (size_t i = 0; i < N; ++i) {
//         for (size_t j = 0; j < N; ++j) {
//             M(i, j) = matrix[i][j];
//         }
//     }
//     return M;
// }

// Eigen::VectorXd VectortoEigenVector(const std::vector<double>& vector) {
//     size_t N = vector.size();
//     Eigen::VectorXd y(N);
//     for (size_t row = 0; row < N; ++row) {
//         y(row) = vector[row];
//     }
//     return V;
// }

// Eigen::MatrixXd compute_V(
//     const std::vector<std::vector<double>>& kinship,
//     double sigma_g_sq,
//     double sigma_e_sq) {

//     Eigen::MatrixXd K = toEigenMatrix(kinship);
//     size_t N = K.rows();
//     Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);
//     return sigma_g_sq * K + sigma_e_sq * I;
// }

// // Compute beta_hat = (X^T V^{-1} X)^{-1} X^T V^{-1} y
// Eigen::VectorXd compute_beta(
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXd& X,
//     const Eigen::VectorXd& y) {

//     // Compute V inverse (use Cholesky for efficiency and stability)
//     Eigen::LLT<Eigen::MatrixXd> lltOfV(V);
//     if(lltOfV.info() != Eigen::Success) {
//         throw std::runtime_error("V matrix decomposition failed");
//     }
//     Eigen::MatrixXd V_inv = lltOfV.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));

//     Eigen::MatrixXd Xt_Vinv = X.transpose() * V_inv;
//     Eigen::MatrixXd Xt_Vinv_X = Xt_Vinv * X;

//     Eigen::VectorXd Xt_Vinv_y = Xt_Vinv * y;

//     // Solve for beta_hat
//     Eigen::VectorXd beta_hat = Xt_Vinv_X.ldlt().solve(Xt_Vinv_y);

//     return beta_hat;
// }

// double compute_reml_log_likelihood(
//     const Eigen::MatrixXd& V,
//     const Eigen::MatrixXd& X,
//     const Eigen::VectorXd& y) {

//     const int N = V.rows();
//     const int p = X.cols();

//     // Cholesky decomposition of V
//     Eigen::LLT<Eigen::MatrixXd> lltOfV(V);
//     if (lltOfV.info() != Eigen::Success) {
//         throw std::runtime_error("V matrix is not positive definite");
//     }

//     // Compute V inverse using Cholesky solve
//     Eigen::MatrixXd V_inv = lltOfV.solve(Eigen::MatrixXd::Identity(N, N));

//     // Compute beta_hat
//     Eigen::MatrixXd Xt_Vinv = X.transpose() * V_inv;
//     Eigen::MatrixXd Xt_Vinv_X = Xt_Vinv * X;
//     Eigen::VectorXd Xt_Vinv_y = Xt_Vinv * y;

//     // Solve for beta_hat
//     Eigen::VectorXd beta_hat = Xt_Vinv_X.ldlt().solve(Xt_Vinv_y);

//     // Compute residuals: y - X beta_hat
//     Eigen::VectorXd resid = y - X * beta_hat;

//     // Compute log determinant of V using Cholesky
//     // log|V| = 2 * sum of log diagonal elements of L, where V = L L^T
//     const auto& L = lltOfV.matrixL();
//     double log_det_V = 0.0;
//     for (int i = 0; i < N; ++i) {
//         log_det_V += std::log(L(i, i));
//     }
//     log_det_V *= 2.0;

//     // Compute log determinant of Xt V^{-1} X
//     Eigen::LDLT<Eigen::MatrixXd> ldlt_XtVinvX(Xt_Vinv_X);
//     if (ldlt_XtVinvX.info() != Eigen::Success) {
//         throw std::runtime_error("Xt_Vinv_X matrix decomposition failed");
//     }
//     double log_det_XtVinvX = 0.0;
//     Eigen::MatrixXd D = ldlt_XtVinvX.vectorD().asDiagonal();
//     for (int i = 0; i < p; ++i) {
//         double val = ldlt_XtVinvX.vectorD()[i];
//         if (val <= 0)
//             throw std::runtime_error("Non-positive diagonal element in Xt_Vinv_X decomposition");
//         log_det_XtVinvX += std::log(val);
//     }

//     // Compute quadratic form resid^T V^{-1} resid
//     double quad_form = resid.transpose() * V_inv * resid;

//     // Compute REML log-likelihood
//     double reml = -0.5 * (log_det_V + log_det_XtVinvX + quad_form + (N - p) * std::log(2.0 * M_PI));

//     return reml;
// }

// struct VarianceComponents {
//     double sigma_g_sq;
//     double sigma_e_sq;
// };

// // Simple optimizer loop to maximize REML over variance components
// VarianceComponents optimize_variance_components(
//     const std::vector<std::vector<double>>& kinship,
//     const Eigen::MatrixXd& X,
//     const Eigen::VectorXd& y,
//     double init_sigma_g_sq = 0.5,
//     double init_sigma_e_sq = 0.5,
//     int max_iter = 100,
//     double tol = 1e-5) {

//     double sigma_g_sq = init_sigma_g_sq;
//     double sigma_e_sq = init_sigma_e_sq;

//     double step = 0.01;  // step size for coordinate ascent
//     double prev_reml = -std::numeric_limits<double>::infinity();

//     for (int iter = 0; iter < max_iter; ++iter) {
//         // --- Optimize sigma_g_sq fixing sigma_e_sq ---
//         double best_sigma_g = sigma_g_sq;
//         double best_reml = prev_reml;

//         // Try small increments and decrements
//         for (double candidate : {sigma_g_sq - step, sigma_g_sq, sigma_g_sq + step}) {
//             if (candidate <= 0) continue;

//             Eigen::MatrixXd V = compute_V(kinship, candidate, sigma_e_sq);
//             double reml;
//             try {
//                 reml = compute_reml_log_likelihood(V, X, y);
//             } catch (...) {
//                 continue;
//             }

//             if (reml > best_reml) {
//                 best_reml = reml;
//                 best_sigma_g = candidate;
//             }
//         }
//         sigma_g_sq = best_sigma_g;
//         prev_reml = best_reml;

//         // --- Optimize sigma_e_sq fixing sigma_g_sq ---
//         double best_sigma_e = sigma_e_sq;
//         best_reml = prev_reml;

//         for (double candidate : {sigma_e_sq - step, sigma_e_sq, sigma_e_sq + step}) {
//             if (candidate <= 0) continue;

//             Eigen::MatrixXd V = compute_V(kinship, sigma_g_sq, candidate);
//             double reml;
//             try {
//                 reml = compute_reml_log_likelihood(V, X, y);
//             } catch (...) {
//                 continue;
//             }

//             if (reml > best_reml) {
//                 best_reml = reml;
//                 best_sigma_e = candidate;
//             }
//         }
//         sigma_e_sq = best_sigma_e;

//         // Check convergence
//         if (std::abs(best_reml - prev_reml) < tol) {
//             break;
//         }
//         prev_reml = best_reml;

//         std::cout << "Iter " << iter << ": sigma_g^2=" << sigma_g_sq
//                   << ", sigma_e^2=" << sigma_e_sq
//                   << ", REML=" << best_reml << std::endl;
//     }

//     return {sigma_g_sq, sigma_e_sq};
// }

// // Logistic link functions
// double logistic(double eta) {
//     return 1.0 / (1.0 + std::exp(-eta));
// }

// // Derivative of logistic inverse (mu) wrt eta
// double logistic_derivative(double eta) {
//     double p = logistic(eta);
//     return p * (1 - p);
// }

// // PQL iteration for binary trait GLMM
// void pql_iteration(
//     const std::vector<std::vector<double>>& kinship,
//     const Eigen::MatrixXd& X,
//     const Eigen::VectorXd& y,
//     Eigen::VectorXd& beta,
//     Eigen::VectorXd& u,
//     double& sigma_g_sq,
//     double& sigma_e_sq,
//     int max_iter = 10,
//     double tol = 1e-5) {

//     const int N = y.size();
//     Eigen::VectorXd eta = X * beta + u; // linear predictor
//     Eigen::VectorXd mu(N);
//     Eigen::VectorXd W_diag(N); // weights
//     Eigen::VectorXd z(N); // working response

//     for (int iter = 0; iter < max_iter; ++iter) {
//         // Step 1: compute mu and weights
//         for (int i = 0; i < N; ++i) {
//             mu[i] = logistic(eta[i]);
//             double dmu_deta = logistic_derivative(eta[i]);
//             // variance for Bernoulli: mu_i * (1 - mu_i)
//             W_diag[i] = dmu_deta * dmu_deta / (mu[i] * (1 - mu[i]) + 1e-6); // avoid div by zero
//             z[i] = eta[i] + (y[i] - mu[i]) / dmu_deta;
//         }

//         // Step 2: transform data by sqrt(W)
//         Eigen::VectorXd sqrt_W = W_diag.array().sqrt();
//         Eigen::MatrixXd X_tilde = X;
//         Eigen::VectorXd z_tilde = z;
//         for (int i = 0; i < N; ++i) {
//             X_tilde.row(i) *= sqrt_W[i];
//             z_tilde[i] *= sqrt_W[i];
//         }

//         // Step 3: compute V matrix with current sigma_g_sq and sigma_e_sq
//         Eigen::MatrixXd V = compute_V(kinship, sigma_g_sq, sigma_e_sq);

//         // Apply weights: V_tilde = W^{1/2} V W^{1/2}
//         // This is approximate, often we treat weights as part of residual variance
//         // For simplicity, just multiply rows and columns by sqrt_W
//         for (int i = 0; i < N; ++i) {
//             for (int j = 0; j < N; ++j) {
//                 V(i, j) *= sqrt_W[i] * sqrt_W[j];
//             }
//         }

//         // Step 4: compute beta update by solving weighted LMM: z_tilde = X_tilde * beta + u + error
//         // For PQL, often random effects are absorbed in V.
//         // Here we solve beta_hat = (X^T V^{-1} X)^{-1} X^T V^{-1} z_tilde
//         try {
//             beta = compute_beta(V, X_tilde, z_tilde);
//         } catch (const std::exception& e) {
//             std::cerr << "Beta estimation failed: " << e.what() << std::endl;
//             return;
//         }

//         // Step 5: update eta and check convergence
//         Eigen::VectorXd eta_new = X * beta + u; // no update for u here, needs more work

//         if ((eta_new - eta).norm() < tol) {
//             std::cout << "PQL converged at iteration " << iter << std::endl;
//             break;
//         }
//         eta = eta_new;

//         // Variance components update can be added here via REML on working LMM
//         // sigma_g_sq, sigma_e_sq = optimize_variance_components(...) using (X_tilde, z_tilde)
//     }
// }

// void lmm_binary(
//     const std::vector<std::vector<double>>& df,              // N x P (paths)
//     const std::vector<bool>& phenotype_binary,               // N
//     const std::vector<std::vector<double>>& kinship,         // N x N
//     const std::vector<std::vector<double>>& covariates,      // N x C
//     std::string& p_value_str, std::string& beta_str,
//     std::string& se_str, std::string& r2_str) {

//     const int N = phenotype_binary.size();
//     const int num_paths = df[0].size();
//     const int num_cov = covariates[0].size();

//     // Convert phenotype std::vector<bool> to Eigen::VectorXd (0/1)
//     Eigen::VectorXd y(N);
//     for (int i = 0; i < N; ++i) y[i] = phenotype_binary[i] ? 1.0 : 0.0;

//     // Convert covariates to Eigen matrix
//     Eigen::MatrixXd cov_mat(N, num_cov);
//     for (int i = 0; i < N; ++i)
//         for (int j = 0; j < num_cov; ++j)
//             cov_mat(i, j) = covariates[i][j];

//     // Build design matrix X = [covariates | df paths]
//     Eigen::MatrixXd X(N, num_cov + num_paths);
//     X.block(0, 0, N, num_cov) = cov_mat;
//     for (int i = 0; i < N; ++i)
//         for (int j = 0; j < num_paths; ++j)
//             X(i, num_cov + j) = df[i][j];

//     // Kinship matrix Eigen conversion
//     Eigen::MatrixXd kinship_mat = MatrixtoEigenMatrix(kinship);

//     // Step 1: Initialize variance components
//     double sigma_g_sq = 0.5, sigma_e_sq = 0.5;

//     // Step 2: Null model optimization (covariates only)
//     VarianceComponents varcomp = optimize_variance_components(
//         kinship, cov_mat, y, sigma_g_sq, sigma_e_sq);
//     sigma_g_sq = varcomp.sigma_g_sq;
//     sigma_e_sq = varcomp.sigma_e_sq;

//     // Step 3: Initialize beta and random effect vector u (zero)
//     Eigen::VectorXd beta = Eigen::VectorXd::Zero(num_cov + num_paths);
//     Eigen::VectorXd u = Eigen::VectorXd::Zero(N);

//     // Step 4: Run PQL iterations to fit GLMM approx for binary trait
//     pql_iteration(kinship, X, y, beta, u, sigma_g_sq, sigma_e_sq);

//     // Step 5: Compute standard errors of betas from final variance components
//     Eigen::MatrixXd V = compute_V(kinship_mat, sigma_g_sq, sigma_e_sq);
//     Eigen::LLT<Eigen::MatrixXd> lltOfV(V);
//     if (lltOfV.info() != Eigen::Success) {
//         throw std::runtime_error("Failed Cholesky decomposition of V");
//     }

//     Eigen::MatrixXd V_inv = lltOfV.solve(Eigen::MatrixXd::Identity(N, N));
//     Eigen::MatrixXd Xt_Vinv = X.transpose() * V_inv;
//     Eigen::MatrixXd Xt_Vinv_X = Xt_Vinv * X;

//     Eigen::MatrixXd cov_beta = Xt_Vinv_X.ldlt().solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
//     Eigen::VectorXd se = cov_beta.diagonal().array().sqrt();

//     // Step 6: Compute p-values and r2 for paths only (skip covariates)
//     std::stringstream pval_ss, beta_ss, se_ss, r2_ss;

//     // Phenotype variance approx
//     double p = y.mean();
//     double var_y = p * (1 - p);

//     for (int path_idx = 0; path_idx < num_paths; ++path_idx) {
//         double b = beta[num_cov + path_idx];
//         double s = se[num_cov + path_idx];
//         if (s == 0) s = 1e-10; // avoid div by zero

//         double z = b / s;
//         double pval = 2 * (1 - std::erf(std::fabs(z) / std::sqrt(2)));

//         // Variance of path allele count
//         Eigen::VectorXd path_vec(N);
//         for (int i = 0; i < N; ++i)
//             path_vec[i] = df[i][path_idx];
//         double var_g = (path_vec.array() - path_vec.mean()).square().mean();

//         double r2 = (b * b * var_g) / var_y;

//         pval_ss << pval << (path_idx == num_paths - 1 ? "" : "\t");
//         beta_ss << b << (path_idx == num_paths - 1 ? "" : "\t");
//         se_ss << s << (path_idx == num_paths - 1 ? "" : "\t");
//         r2_ss << r2 << (path_idx == num_paths - 1 ? "" : "\t");
//     }

//     p_value_str = pval_ss.str();
//     beta_str = beta_ss.str();
//     se_str = se_ss.str();
//     r2_str = r2_ss.str();
// }

// void lmm_binary(
//     const std::vector<std::vector<double>>& df,              // N x P (paths)
//     const std::vector<bool>& phenotype_binary,               // N
//     const stoat_vcf::KinshipMatrix& kinship,                                              
//     const std::vector<std::vector<double>>& covariates,      // N x C
//     std::string& p_value_str, std::string& beta_str,
//     std::string& se_str, std::string& r2_str) {
// }

// void lmm_quantitative(
//     const std::vector<std::vector<double>>& df,                  
//     const std::vector<double>& phenotype_table,      
//     const stoat_vcf::KinshipMatrix& kinship,                                              
//     const std::vector<std::vector<double>>& covariates,
//     std::string& p_value_str, std::string& beta_str, 
//     std::string& se_str, std::string& r2_str) {
// }
