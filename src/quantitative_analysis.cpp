#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"

using namespace std;

// Linear regression function OLS with intercept
void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {

    size_t num_samples = df.size();
    size_t max_paths = df[0].size();

    Eigen::MatrixXd X(num_samples, max_paths);
    X.setZero(); // Initialize matrix with zeros
    Eigen::VectorXd y(num_samples);

    for (size_t row=0; row < num_samples; ++row) {
        y(row) = quantitative_phenotype[row];
        for (size_t col = 0; col < max_paths; ++col) {
            X(row, col) = static_cast<double>(df[row][col]);
        }
    }

    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;
    
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = max_paths - 1; // degrees of freedom
    int df_res = num_samples - max_paths;
    double mse = rss / df_res;  // Mean Squared Error (MSE)
    
    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // Compute F-statistic
    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    boost::math::fisher_f dist(df_reg, df_res);
    double p_value = boost::math::cdf(boost::math::complement(dist, std::abs(f_stat)));

    // set precision : 4 digit
    r2_str = set_precision(r2);
    beta_str = set_precision(beta.mean());
    se_str = set_precision(se.mean());
    p_value_str = set_precision(p_value);
}

void glm_quantitative(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = covar.empty() ? 0 : covar[0].size();
    size_t num_features = num_variants + num_covariates;

    Eigen::MatrixXd X(num_samples, num_features);
    Eigen::VectorXd y(num_samples);
    
    for (size_t i = 0; i < num_samples; ++i) {
        y(i) = quantitative_phenotype[i];
        size_t col = 0;
        for (size_t j = 0; j < num_variants; ++j) {
            X(i, col++) = df[i][j];
        }
        for (size_t j = 0; j < num_covariates; ++j) {
            X(i, col++) = covar[i][j];
        }
    }
    
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1 - (rss / tss);

    int df_reg = num_features - 1;
    int df_res = num_samples - num_features;
    double mse = rss / df_res;

    Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    double f_stat = (r2 / df_reg) / ((1 - r2) / df_res);
    boost::math::fisher_f dist(df_reg, df_res);
    double p_value = boost::math::cdf(boost::math::complement(dist, std::abs(f_stat)));

    // Output (mean of all variant betas and SEs)
    beta_str = set_precision(beta.segment(1, num_variants).mean());
    se_str = set_precision(se.segment(1, num_variants).mean());
    r2_str = set_precision(r2);
    p_value_str = set_precision(p_value);
}

// Explicit template instantiations
template std::tuple<std::vector<std::vector<double>>, std::vector<double>, size_t>
create_quantitative_table<double>(
    const size_t&,
    const std::vector<std::string>&,
    const std::vector<double>&,
    Matrix&);

template std::tuple<std::vector<std::vector<double>>, std::vector<bool>, size_t>
create_quantitative_table<bool>(
    const size_t&,
    const std::vector<std::string>&,
    const std::vector<bool>&,
    Matrix&);

// Function template definition
template <typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>, size_t> create_quantitative_table(
    const size_t& length_sample,
    const std::vector<std::string>& column_headers,
    const std::vector<T>& phenotype,
    Matrix& matrix) {

    size_t allele_number = 0;
    size_t length_column = column_headers.size();

    std::vector<std::vector<double>> genotypes(length_sample, std::vector<double>(length_column, 0.0));
    std::unordered_set<size_t> index_used;

    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        std::vector<size_t> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, length_sample * 2);

        for (size_t idx : idx_srr_save) {
            size_t srr_idx = idx / 2;
            genotypes[srr_idx][col_idx] += 1.0;
            index_used.insert(srr_idx);
            allele_number++;
        }
    }

    std::vector<std::vector<double>> genotypes_filtered;
    genotypes_filtered.reserve(index_used.size());

    std::vector<T> phenotype_filtered;
    phenotype_filtered.reserve(index_used.size());

    for (size_t i : index_used) {
        double row_sum = std::accumulate(genotypes[i].begin(), genotypes[i].end(), 0.0);

        std::vector<double> normalized_row;
        normalized_row.reserve(length_column);
        for (double allele : genotypes[i]) {
            normalized_row.push_back(allele > 0.0 ? allele / row_sum : 0.0);
        }

        genotypes_filtered.push_back(std::move(normalized_row));
        phenotype_filtered.push_back(phenotype[i]);
    }

    return {genotypes_filtered, phenotype_filtered, allele_number};
}

std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, size_t> create_eqtl_table(
    const size_t& length_sample,
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    size_t allele_number = 0;
    size_t length_column = column_headers.size();

    std::vector<std::vector<double>> genotypes(length_sample, std::vector<double>(length_column, 0.0));
    std::unordered_set<size_t> index_used;

    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        std::vector<size_t> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, length_sample * 2);

        for (size_t idx : idx_srr_save) {
            size_t srr_idx = idx / 2;
            genotypes[srr_idx][col_idx] += 1.0;
            index_used.insert(srr_idx);
            allele_number++;
        }
    }

    std::vector<std::vector<double>> genotypes_filtered;
    genotypes_filtered.reserve(index_used.size());

    for (size_t i : index_used) {
        double row_sum = std::accumulate(genotypes[i].begin(), genotypes[i].end(), 0.0);

        std::vector<double> normalized_row;
        normalized_row.reserve(length_column);
        for (double allele : genotypes[i]) {
            normalized_row.push_back(allele > 0.0 ? allele / row_sum : 0.0);
        }
        genotypes_filtered.push_back(std::move(normalized_row));
    }

    return {genotypes_filtered, index_used, allele_number};
}