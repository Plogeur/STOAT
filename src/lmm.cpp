#include "lmm.hpp"

// Function to calculate the Normal CDF (for p-values)
double normal_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2));
}

// Updated LMM function accepting new data structures
std::tuple<string, string, string, string> lmm_quantitative(
    const std::unordered_map<std::string, std::vector<int>>& df,                  
    const std::unordered_map<std::string, double>& phenotype_table,      
    const KinshipMatrix& kinship,                                              
    const std::unordered_map<std::string, std::vector<double>>& covariates) {

    const int N = kinship.ids.size();  // Number of samples

    // Prepare phenotype vector
    Eigen::VectorXd phenotype(N);
    for (int i = 0; i < N; ++i) {
        phenotype(i) = phenotype_table.at(kinship.ids[i]);
    }

    // Prepare kinship matrix
    Eigen::MatrixXd kinship_matrix(N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            kinship_matrix(i, j) = kinship.matrix[i][j];
        }
    }

    // Step 2: Prepare genotype matrix (X matrix)
    Eigen::MatrixXd X(N, 1);  // For simplicity, using only 1 SNP here
    for (int i = 0; i < N; ++i) {
        const std::string& sample_id = kinship.ids[i];
        const auto& allele_counts = df.at(sample_id);
        int genotype = std::accumulate(allele_counts.begin(), allele_counts.end(), 0);
        X(i, 0) = static_cast<double>(genotype);  // SNP dosage: 0, 1, 2
    }

    // Prepare covariates matrix (additional columns)
    Eigen::MatrixXd cov(N, covariates.begin()->second.size());
    for (int i = 0; i < N; ++i) {
        const std::string& sample_id = kinship.ids[i];
        const auto& covariate_vals = covariates.at(sample_id);
        for (int j = 0; j < covariate_vals.size(); ++j) {
            cov(i, j) = covariate_vals[j];
        }
    }

    // Combine SNP genotype (X) with covariates matrix
    Eigen::MatrixXd X_full(N, cov.cols() + 1);
    X_full << X, cov;  // First column is the SNP, rest are covariates

    // Step 3: Perform LMM calculations
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(kinship_matrix);
    Eigen::MatrixXd U = eig.eigenvectors();
    Eigen::VectorXd S = eig.eigenvalues();

    Eigen::VectorXd y_star = U.transpose() * phenotype;
    Eigen::MatrixXd X_star = U.transpose() * X_full;

    // NOTE: This delta should ideally be estimated (REML) in real models
    double delta = 1.0;

    // Step 4: Calculate inverse variance
    Eigen::VectorXd V_diag = S.array() + delta;
    Eigen::VectorXd V_inv_diag = V_diag.cwiseInverse();

    // Step 5: GLS (Generalized Least Squares) estimation
    Eigen::MatrixXd XVX = X_star.transpose() * V_inv_diag.asDiagonal() * X_star;
    Eigen::VectorXd XVy = X_star.transpose() * V_inv_diag.asDiagonal() * y_star;

    Eigen::VectorXd beta_hat = XVX.ldlt().solve(XVy);

    // Calculate standard error and t-statistic for SNP effect (first beta)
    Eigen::MatrixXd XVX_inv = XVX.inverse();
    double beta = beta_hat(0);
    double se = std::sqrt(XVX_inv(0, 0));
    double t_stat = beta / se;
    double p_val = 2.0 * (1.0 - normal_cdf(std::fabs(t_stat)));

    string beta_str = set_precision(beta);
    string se_str = set_precision(se);
    string t_stat_str = set_precision(t_stat);
    string p_val_str = set_precision(p_val);

    return {t_stat_str, beta_str, se_str, p_val_str};
}

