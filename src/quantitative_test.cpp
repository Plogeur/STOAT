#include "quatitative_test.hpp"
#include "snarl_analyser.hpp"
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
    size_t num_features = df[0].size();

    // Add intercept: X with one additional column for intercept
    Eigen::MatrixXd X(num_samples, num_features + 1);
    X.col(0) = Eigen::VectorXd::Ones(num_samples);  // Intercept column
    for (size_t row = 0; row < num_samples; ++row) {
        for (size_t col = 0; col < num_features; ++col) {
            X(row, col + 1) = df[row][col];
        }
    }

    // Response vector
    Eigen::VectorXd y(num_samples);
    for (size_t row = 0; row < num_samples; ++row) {
        y(row) = quantitative_phenotype[row];
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
    for (int i = 1; i < num_features+1; ++i) { // i = 1 avoid const p-value
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

// Linear regression function OLS with intercept + covariate
void glm_quantitative(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str) {

    size_t num_samples = df.size();
    size_t num_variants = df[0].size();
    size_t num_covariates = covar[0].size();
    size_t num_features = num_variants + num_covariates;

    Eigen::MatrixXd X(num_samples, num_features + 1);
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
    
    Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // R²
    double rss = residuals.squaredNorm();
    double tss = (y.array() - y.mean()).matrix().squaredNorm();
    double r2 = 1.0 - (rss / tss);

    int df_res = num_samples - X.cols();    // residual degrees of freedom
    double mse = rss / df_res;

    // Standard errors
    Eigen::LDLT<Eigen::MatrixXd> ldlt(X.transpose() * X);
    Eigen::MatrixXd cov_matrix = ldlt.solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
    Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt();

    // Standard errors
    // Eigen::MatrixXd cov_matrix = (X.transpose() * X).inverse();
    // Eigen::VectorXd se = (cov_matrix.diagonal() * mse).array().sqrt().matrix();

    // t-statistics
    Eigen::VectorXd t_stats = beta.array() / se.array();
    boost::math::students_t t_dist(df_res);
    std::vector<double> p_values;
    for (int i = 1; i < num_variants+1; ++i) {
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

// Explicit template instantiations
template std::tuple<std::vector<std::vector<double>>, std::vector<double>, size_t, std::vector<size_t>>
create_quantitative_table<double>(
    const size_t&,
    const std::vector<Path_traversal_t>&,
    const std::vector<double>&,
    const EdgeBySampleMatrix&);

template std::tuple<std::vector<std::vector<double>>, std::vector<bool>, size_t, std::vector<size_t>>
create_quantitative_table<bool>(
    const size_t&,
    const std::vector<Path_traversal_t>&,
    const std::vector<bool>&,
    const EdgeBySampleMatrix&);

std::tuple<std::vector<std::vector<double>>, size_t, std::unordered_set<size_t>, bool, std::vector<size_t>>
    process_table_quantitative(
        const size_t& number_samples,
        const std::vector<Path_traversal_t>& column_headers,
        const EdgeBySampleMatrix& matrix) {

    size_t allele_number = 0;
    size_t length_column = column_headers.size();

    std::vector<size_t> allele_paths(length_column, 0);

    std::vector<std::vector<double>> genotypes(number_samples);
    for (auto& row : genotypes)
        row.reserve(length_column);

    std::vector<size_t> kept_columns; // Indices of valid columns
    std::unordered_set<size_t> index_used;

    // Loop over all columns
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const Path_traversal_t& path_snarl = column_headers[col_idx];
        std::vector<Edge_t> list_edge_path = decompose_path_to_edges(path_snarl);

        //Get the indices of all samples that take this path
        std::vector<size_t> idx_srr_save = identify_path(list_edge_path, matrix, number_samples * 2);

        if (idx_srr_save.empty())
            continue; // Skip if column is empty

        kept_columns.push_back(col_idx); // Valid column

        // Ensure rows have space for new column
        for (size_t i = 0; i < number_samples; ++i) {
            if (genotypes[i].size() < kept_columns.size())
                genotypes[i].resize(kept_columns.size(), 0.0);
        }

        size_t numb_all = idx_srr_save.size();
        allele_number += numb_all;
        allele_paths[col_idx] = numb_all;

        // Fill genotype matrix
        for (size_t idx : idx_srr_save) {
            size_t srr_idx = idx / 2;
            genotypes[srr_idx][kept_columns.size() - 1] += 1.0;
            index_used.insert(srr_idx);
        }
    }

    // Trim last column if needed
    bool drop_last_col = (kept_columns.size() > 1); 

    return {genotypes, allele_number, index_used, drop_last_col, allele_paths};   
}

// Function template definition
template<typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>, size_t, std::vector<size_t>> create_quantitative_table(
    const size_t& number_samples,
    const std::vector<Path_traversal_t>& column_headers,
    const std::vector<T>& phenotype,
    const EdgeBySampleMatrix& matrix) {

    const auto& [genotypes, allele_number, index_used, drop_last_col, allele_paths] = 
    process_table_quantitative(number_samples, column_headers, matrix);

    std::vector<std::vector<double>> genotypes_filtered;
    genotypes_filtered.reserve(index_used.size());

    std::vector<T> phenotype_filtered;
    phenotype_filtered.reserve(index_used.size());

    for (size_t i : index_used) {
        const auto& row = genotypes[i];
        double row_sum = std::accumulate(row.begin(), row.end(), 0.0);

        std::vector<double> normalized_row;
        size_t max_col = drop_last_col ? row.size() - 1 : row.size();
        normalized_row.reserve(max_col);

        for (size_t j = 0; j < max_col; ++j) {
            normalized_row.push_back(row[j] > 0.0 ? row[j] / row_sum : 0.0);
        }

        genotypes_filtered.push_back(std::move(normalized_row));
        phenotype_filtered.push_back(phenotype[i]);
    }

    return {genotypes_filtered, phenotype_filtered, allele_number, allele_paths};
}

std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, size_t, std::vector<size_t>> create_eqtl_table(
    const size_t& number_samples,
    const std::vector<Path_traversal_t>& column_headers,
    const EdgeBySampleMatrix& matrix) {

    const auto& [genotypes, allele_number, index_used, drop_last_col, allele_paths] = 
    process_table_quantitative(number_samples, column_headers, matrix);

    std::vector<std::vector<double>> genotypes_filtered;
    genotypes_filtered.reserve(index_used.size());

    for (size_t i : index_used) {
        const auto& row = genotypes[i];
        double row_sum = std::accumulate(row.begin(), row.end(), 0.0);

        std::vector<double> normalized_row;
        size_t max_col = drop_last_col ? row.size() - 1 : row.size();
        normalized_row.reserve(max_col);

        for (size_t j = 0; j < max_col; ++j) {
            normalized_row.push_back(row[j] > 0.0 ? row[j] / row_sum : 0.0);
        }

        genotypes_filtered.push_back(std::move(normalized_row));
    }

    return {genotypes_filtered, index_used, allele_number, allele_paths};
}
