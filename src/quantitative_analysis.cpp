#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

using namespace std;

// Linear regression function OLS
std::tuple<std::string, std::string, std::string, std::string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {

    size_t num_samples = df.size();
    size_t max_paths = 0;
    
    // Determine the maximum number of predictors (columns)
    for (const auto& [_, paths] : df) {
        max_paths = std::max(max_paths, paths.size());
    }
    
    Eigen::MatrixXd X(num_samples, max_paths);
    X.setZero(); // Initialize matrix with zeros
    Eigen::VectorXd y(num_samples);
    
    int row = 0;
    for (const auto& [sample, paths] : df) {
        y(row) = quantitative_phenotype.at(sample);
        for (size_t col = 0; col < paths.size(); ++col) {
            X(row, col) = paths[col];
        }
        ++row;
    }
    
    // Perform OLS: β = (X^T X)^(-1) X^T y
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    
    // Compute residuals
    Eigen::VectorXd residuals = y - X * beta;
    double ss_res = residuals.squaredNorm();
    
    // Compute total sum of squares
    double y_mean = y.mean();
    double ss_tot = (y.array() - y_mean).square().sum();
    
    // Compute R-squared
    double r_squared = 1 - (ss_res / ss_tot);
    
    // Compute F-statistic
    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    
    double p_value = 1.0f;
    if (f_stat > 0) {
        boost::math::fisher_f dist(df_reg, df_res);
        p_value = boost::math::cdf(boost::math::complement(dist, f_stat));
    }

    // set precision : 4 digit
    string r2_str = set_precision(r2);
    string beta_mean_str = set_precision(beta.mean());
    string se_mean_str = set_precision(se.mean());
    string p_value_str = set_precision(global_p_value);

    return {r2_str, beta_mean_str, se_mean_str, p_value_str};
}

// Function to fit a Linear Mixed Model (LMM) for GWAS and return global p-value, mean beta, and mean SE
std::tuple<std::string, std::string, std::string> LMM_quantitative(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covariate) {

    // Convert phenotype data into an Eigen vector
    Eigen::VectorXd y(quantitative_phenotype.size());
    Eigen::MatrixXd covariate_matrix(quantitative_phenotype.size(), covariate.begin()->second.size());
    int index = 0;
    std::vector<std::string> sample_ids;

    for (const auto& pair : quantitative_phenotype) {
        y(index) = pair.second;
        sample_ids.push_back(pair.first);
        for (int j = 0; j < covariate.at(pair.first).size(); j++) {
            covariate_matrix(index, j) = covariate.at(pair.first)[j];
        }
        index++;
    }

    // Construct the genotype matrix X
    Eigen::MatrixXd X(sample_ids.size(), df.size());
    int col = 0;
    for (const auto& snp : df) {
        for (int row = 0; row < sample_ids.size(); row++) {
            X(row, col) = snp.second[row];
        }
        col++;
    }

    // Combine genotype matrix with covariates
    Eigen::MatrixXd design_matrix(sample_ids.size(), X.cols() + covariate_matrix.cols());
    design_matrix << X, covariate_matrix;

    // Estimate beta using OLS (beta = (X'X)^-1 X'Y)
    Eigen::MatrixXd XtX = design_matrix.transpose() * design_matrix;
    Eigen::VectorXd Xty = design_matrix.transpose() * y;
    Eigen::VectorXd beta = XtX.ldlt().solve(Xty);

    // Compute residuals
    Eigen::VectorXd residuals = y - design_matrix * beta;
    double residual_variance = (residuals.transpose() * residuals).value() / residuals.size();

    // Compute the covariance matrix of the estimated coefficients
    Eigen::MatrixXd XtX_inv = XtX.ldlt().solve(Eigen::MatrixXd::Identity(XtX.rows(), XtX.cols()));
    Eigen::MatrixXd beta_cov = residual_variance * XtX_inv;

    // Compute standard errors (diagonal of the covariance matrix)
    Eigen::VectorXd standard_errors = beta_cov.diagonal().array().sqrt();

    // Compute F-statistic for the global test (based on the model's explained variance)
    Eigen::VectorXd y_hat = design_matrix * beta; // predicted values
    Eigen::VectorXd residuals_full = y - y_hat; // residuals
    double ss_total = (y - y.mean()).transpose() * (y - y.mean());
    double ss_residual = residuals_full.transpose() * residuals_full;
    double ss_model = ss_total - ss_residual;

    int degrees_of_freedom_model = X.cols() + covariate_matrix.cols();
    int degrees_of_freedom_residual = y.size() - degrees_of_freedom_model;

    double f_statistic = (ss_model / degrees_of_freedom_model) / (ss_residual / degrees_of_freedom_residual);
    
    // Compute p-value from the F-statistic (using the F-distribution)
    double global_p_value = 1 - std::exp(-f_statistic / 2); // Simplified approximation for the p-value

    // Convert the global p-value to a string
    std::ostringstream p_value_stream;
    p_value_stream << global_p_value;
    std::string p_value_str = p_value_stream.str();

    // Compute mean beta and mean standard error
    double mean_beta = beta.mean();
    double mean_se = standard_errors.mean();

    // Convert mean beta and mean SE to strings
    std::ostringstream beta_stream;
    beta_stream << mean_beta;
    std::string beta_str = beta_stream.str();

    std::ostringstream se_stream;
    se_stream << mean_se;
    std::string se_str = se_stream.str();

    // Return the tuple with p-value, mean beta, and mean SE as strings
    return std::make_tuple(p_value_str, beta_str, se_str);
}

// Function to create the quantitative table
std::pair<std::unordered_map<std::string, std::vector<int>>, size_t> create_quantitative_table(
    const std::vector<std::string>& list_samples, 
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    // Retrieve row headers dictionary
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();
    size_t allele_number = 0;
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
            allele_number++;
        }
    }

    std::unordered_map<std::string, std::vector<int>> df;
    for (size_t i = 0; i < list_samples.size(); ++i) {
        df[list_samples[i]] = genotypes[i];
    }
    
    return {df, allele_number};
}
