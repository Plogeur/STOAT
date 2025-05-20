#include "lmm.hpp"

Eigen::MatrixXd MatrixtoEigenMatrix(const std::vector<std::vector<double>>& matrix) {
    size_t N = matrix.size();
    Eigen::MatrixXd M(N, N);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            M(i, j) = matrix[i][j];
        }
    }
    return M;
}

Eigen::VectorXd VectortoEigenVector(const std::vector<double>& vector) {
    size_t N = vector.size();
    Eigen::VectorXd y(N);
    for (size_t row = 0; row < N; ++row) {
        y(row) = vector[row];
    }
    return V;
}

Eigen::MatrixXd compute_V(
    const std::vector<std::vector<double>>& kinship,
    double sigma_g_sq,
    double sigma_e_sq) {

    Eigen::MatrixXd K = toEigenMatrix(kinship);
    size_t N = K.rows();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);
    return sigma_g_sq * K + sigma_e_sq * I;
}

// Compute beta_hat = (X^T V^{-1} X)^{-1} X^T V^{-1} y
Eigen::VectorXd compute_beta(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y) {

    // Compute V inverse (use Cholesky for efficiency and stability)
    Eigen::LLT<Eigen::MatrixXd> lltOfV(V);
    if(lltOfV.info() != Eigen::Success) {
        throw std::runtime_error("V matrix decomposition failed");
    }
    Eigen::MatrixXd V_inv = lltOfV.solve(Eigen::MatrixXd::Identity(V.rows(), V.cols()));

    Eigen::MatrixXd Xt_Vinv = X.transpose() * V_inv;
    Eigen::MatrixXd Xt_Vinv_X = Xt_Vinv * X;

    Eigen::VectorXd Xt_Vinv_y = Xt_Vinv * y;

    // Solve for beta_hat
    Eigen::VectorXd beta_hat = Xt_Vinv_X.ldlt().solve(Xt_Vinv_y);

    return beta_hat;
}

double compute_reml_log_likelihood(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y) {

    const int N = V.rows();
    const int p = X.cols();

    // Cholesky decomposition of V
    Eigen::LLT<Eigen::MatrixXd> lltOfV(V);
    if (lltOfV.info() != Eigen::Success) {
        throw std::runtime_error("V matrix is not positive definite");
    }

    // Compute V inverse using Cholesky solve
    Eigen::MatrixXd V_inv = lltOfV.solve(Eigen::MatrixXd::Identity(N, N));

    // Compute beta_hat
    Eigen::MatrixXd Xt_Vinv = X.transpose() * V_inv;
    Eigen::MatrixXd Xt_Vinv_X = Xt_Vinv * X;
    Eigen::VectorXd Xt_Vinv_y = Xt_Vinv * y;

    // Solve for beta_hat
    Eigen::VectorXd beta_hat = Xt_Vinv_X.ldlt().solve(Xt_Vinv_y);

    // Compute residuals: y - X beta_hat
    Eigen::VectorXd resid = y - X * beta_hat;

    // Compute log determinant of V using Cholesky
    // log|V| = 2 * sum of log diagonal elements of L, where V = L L^T
    const auto& L = lltOfV.matrixL();
    double log_det_V = 0.0;
    for (int i = 0; i < N; ++i) {
        log_det_V += std::log(L(i, i));
    }
    log_det_V *= 2.0;

    // Compute log determinant of Xt V^{-1} X
    Eigen::LDLT<Eigen::MatrixXd> ldlt_XtVinvX(Xt_Vinv_X);
    if (ldlt_XtVinvX.info() != Eigen::Success) {
        throw std::runtime_error("Xt_Vinv_X matrix decomposition failed");
    }
    double log_det_XtVinvX = 0.0;
    Eigen::MatrixXd D = ldlt_XtVinvX.vectorD().asDiagonal();
    for (int i = 0; i < p; ++i) {
        double val = ldlt_XtVinvX.vectorD()[i];
        if (val <= 0)
            throw std::runtime_error("Non-positive diagonal element in Xt_Vinv_X decomposition");
        log_det_XtVinvX += std::log(val);
    }

    // Compute quadratic form resid^T V^{-1} resid
    double quad_form = resid.transpose() * V_inv * resid;

    // Compute REML log-likelihood
    double reml = -0.5 * (log_det_V + log_det_XtVinvX + quad_form + (N - p) * std::log(2.0 * M_PI));

    return reml;
}

struct VarianceComponents {
    double sigma_g_sq;
    double sigma_e_sq;
};

// Simple optimizer loop to maximize REML over variance components
VarianceComponents optimize_variance_components(
    const std::vector<std::vector<double>>& kinship,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    double init_sigma_g_sq = 0.5,
    double init_sigma_e_sq = 0.5,
    int max_iter = 100,
    double tol = 1e-5) {

    double sigma_g_sq = init_sigma_g_sq;
    double sigma_e_sq = init_sigma_e_sq;

    double step = 0.01;  // step size for coordinate ascent
    double prev_reml = -std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < max_iter; ++iter) {
        // --- Optimize sigma_g_sq fixing sigma_e_sq ---
        double best_sigma_g = sigma_g_sq;
        double best_reml = prev_reml;

        // Try small increments and decrements
        for (double candidate : {sigma_g_sq - step, sigma_g_sq, sigma_g_sq + step}) {
            if (candidate <= 0) continue;

            Eigen::MatrixXd V = compute_V(kinship, candidate, sigma_e_sq);
            double reml;
            try {
                reml = compute_reml_log_likelihood(V, X, y);
            } catch (...) {
                continue;
            }

            if (reml > best_reml) {
                best_reml = reml;
                best_sigma_g = candidate;
            }
        }
        sigma_g_sq = best_sigma_g;
        prev_reml = best_reml;

        // --- Optimize sigma_e_sq fixing sigma_g_sq ---
        double best_sigma_e = sigma_e_sq;
        best_reml = prev_reml;

        for (double candidate : {sigma_e_sq - step, sigma_e_sq, sigma_e_sq + step}) {
            if (candidate <= 0) continue;

            Eigen::MatrixXd V = compute_V(kinship, sigma_g_sq, candidate);
            double reml;
            try {
                reml = compute_reml_log_likelihood(V, X, y);
            } catch (...) {
                continue;
            }

            if (reml > best_reml) {
                best_reml = reml;
                best_sigma_e = candidate;
            }
        }
        sigma_e_sq = best_sigma_e;

        // Check convergence
        if (std::abs(best_reml - prev_reml) < tol) {
            break;
        }
        prev_reml = best_reml;

        std::cout << "Iter " << iter << ": sigma_g^2=" << sigma_g_sq
                  << ", sigma_e^2=" << sigma_e_sq
                  << ", REML=" << best_reml << std::endl;
    }

    return {sigma_g_sq, sigma_e_sq};
}

// Logistic link functions
double logistic(double eta) {
    return 1.0 / (1.0 + std::exp(-eta));
}

// Derivative of logistic inverse (mu) wrt eta
double logistic_derivative(double eta) {
    double p = logistic(eta);
    return p * (1 - p);
}

// PQL iteration for binary trait GLMM
void pql_iteration(
    const std::vector<std::vector<double>>& kinship,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    Eigen::VectorXd& beta,
    Eigen::VectorXd& u,
    double& sigma_g_sq,
    double& sigma_e_sq,
    int max_iter = 10,
    double tol = 1e-5) {

    const int N = y.size();
    Eigen::VectorXd eta = X * beta + u; // linear predictor
    Eigen::VectorXd mu(N);
    Eigen::VectorXd W_diag(N); // weights
    Eigen::VectorXd z(N); // working response

    for (int iter = 0; iter < max_iter; ++iter) {
        // Step 1: compute mu and weights
        for (int i = 0; i < N; ++i) {
            mu[i] = logistic(eta[i]);
            double dmu_deta = logistic_derivative(eta[i]);
            // variance for Bernoulli: mu_i * (1 - mu_i)
            W_diag[i] = dmu_deta * dmu_deta / (mu[i] * (1 - mu[i]) + 1e-6); // avoid div by zero
            z[i] = eta[i] + (y[i] - mu[i]) / dmu_deta;
        }

        // Step 2: transform data by sqrt(W)
        Eigen::VectorXd sqrt_W = W_diag.array().sqrt();
        Eigen::MatrixXd X_tilde = X;
        Eigen::VectorXd z_tilde = z;
        for (int i = 0; i < N; ++i) {
            X_tilde.row(i) *= sqrt_W[i];
            z_tilde[i] *= sqrt_W[i];
        }

        // Step 3: compute V matrix with current sigma_g_sq and sigma_e_sq
        Eigen::MatrixXd V = compute_V(kinship, sigma_g_sq, sigma_e_sq);

        // Apply weights: V_tilde = W^{1/2} V W^{1/2}
        // This is approximate, often we treat weights as part of residual variance
        // For simplicity, just multiply rows and columns by sqrt_W
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                V(i, j) *= sqrt_W[i] * sqrt_W[j];
            }
        }

        // Step 4: compute beta update by solving weighted LMM: z_tilde = X_tilde * beta + u + error
        // For PQL, often random effects are absorbed in V.
        // Here we solve beta_hat = (X^T V^{-1} X)^{-1} X^T V^{-1} z_tilde
        try {
            beta = compute_beta(V, X_tilde, z_tilde);
        } catch (const std::exception& e) {
            std::cerr << "Beta estimation failed: " << e.what() << std::endl;
            return;
        }

        // Step 5: update eta and check convergence
        Eigen::VectorXd eta_new = X * beta + u; // no update for u here, needs more work

        if ((eta_new - eta).norm() < tol) {
            std::cout << "PQL converged at iteration " << iter << std::endl;
            break;
        }
        eta = eta_new;

        // Variance components update can be added here via REML on working LMM
        // sigma_g_sq, sigma_e_sq = optimize_variance_components(...) using (X_tilde, z_tilde)
    }
}

void lmm_binary(
    const std::vector<std::vector<double>>& df,              // N x P (paths)
    const std::vector<bool>& phenotype_binary,               // N
    const std::vector<std::vector<double>>& kinship,         // N x N
    const std::vector<std::vector<double>>& covariates,      // N x C
    std::string& p_value_str, std::string& beta_str,
    std::string& se_str, std::string& r2_str) {

    const int N = phenotype_binary.size();
    const int num_paths = df[0].size();
    const int num_cov = covariates[0].size();

    // Convert phenotype vector<bool> to Eigen::VectorXd (0/1)
    Eigen::VectorXd y(N);
    for (int i = 0; i < N; ++i) y[i] = phenotype_binary[i] ? 1.0 : 0.0;

    // Convert covariates to Eigen matrix
    Eigen::MatrixXd cov_mat(N, num_cov);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < num_cov; ++j)
            cov_mat(i, j) = covariates[i][j];

    // Build design matrix X = [covariates | df paths]
    Eigen::MatrixXd X(N, num_cov + num_paths);
    X.block(0, 0, N, num_cov) = cov_mat;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < num_paths; ++j)
            X(i, num_cov + j) = df[i][j];

    // Kinship matrix Eigen conversion
    Eigen::MatrixXd kinship_mat = MatrixtoEigenMatrix(kinship);

    // Step 1: Initialize variance components
    double sigma_g_sq = 0.5, sigma_e_sq = 0.5;

    // Step 2: Null model optimization (covariates only)
    VarianceComponents varcomp = optimize_variance_components(
        kinship, cov_mat, y, sigma_g_sq, sigma_e_sq);
    sigma_g_sq = varcomp.sigma_g_sq;
    sigma_e_sq = varcomp.sigma_e_sq;

    // Step 3: Initialize beta and random effect vector u (zero)
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(num_cov + num_paths);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(N);

    // Step 4: Run PQL iterations to fit GLMM approx for binary trait
    pql_iteration(kinship, X, y, beta, u, sigma_g_sq, sigma_e_sq);

    // Step 5: Compute standard errors of betas from final variance components
    Eigen::MatrixXd V = compute_V(kinship_mat, sigma_g_sq, sigma_e_sq);
    Eigen::LLT<Eigen::MatrixXd> lltOfV(V);
    if (lltOfV.info() != Eigen::Success) {
        throw std::runtime_error("Failed Cholesky decomposition of V");
    }

    Eigen::MatrixXd V_inv = lltOfV.solve(Eigen::MatrixXd::Identity(N, N));
    Eigen::MatrixXd Xt_Vinv = X.transpose() * V_inv;
    Eigen::MatrixXd Xt_Vinv_X = Xt_Vinv * X;

    Eigen::MatrixXd cov_beta = Xt_Vinv_X.ldlt().solve(Eigen::MatrixXd::Identity(Xt_Vinv_X.rows(), Xt_Vinv_X.cols()));
    Eigen::VectorXd se = cov_beta.diagonal().array().sqrt();

    // Step 6: Compute p-values and r2 for paths only (skip covariates)
    std::stringstream pval_ss, beta_ss, se_ss, r2_ss;

    // Phenotype variance approx
    double p = y.mean();
    double var_y = p * (1 - p);

    for (int path_idx = 0; path_idx < num_paths; ++path_idx) {
        double b = beta[num_cov + path_idx];
        double s = se[num_cov + path_idx];
        if (s == 0) s = 1e-10; // avoid div by zero

        double z = b / s;
        double pval = 2 * (1 - std::erf(std::fabs(z) / std::sqrt(2)));

        // Variance of path allele count
        Eigen::VectorXd path_vec(N);
        for (int i = 0; i < N; ++i)
            path_vec[i] = df[i][path_idx];
        double var_g = (path_vec.array() - path_vec.mean()).square().mean();

        double r2 = (b * b * var_g) / var_y;

        pval_ss << pval << (path_idx == num_paths - 1 ? "" : "\t");
        beta_ss << b << (path_idx == num_paths - 1 ? "" : "\t");
        se_ss << s << (path_idx == num_paths - 1 ? "" : "\t");
        r2_ss << r2 << (path_idx == num_paths - 1 ? "" : "\t");
    }

    p_value_str = pval_ss.str();
    beta_str = beta_ss.str();
    se_str = se_ss.str();
    r2_str = r2_ss.str();
}

void lmm_quantitative(
    const std::vector<std::vector<double>>& df,                  
    const vector<double>& phenotype_table,      
    const KinshipMatrix& kinship,                                              
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {
}

