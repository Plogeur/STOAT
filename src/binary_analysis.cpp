#include "binary_analysis.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

// ------------------------ LMM LOGISTIC REGRESSION ------------------------

// Function to compute p-values from Wald test
Eigen::VectorXd computePValues(const Eigen::VectorXd& beta, const Eigen::MatrixXd& XtWX) {
    Eigen::VectorXd standard_errors = XtWX.diagonal().array().sqrt().inverse();
    Eigen::VectorXd z_scores = beta.array() / standard_errors.array();
    Eigen::VectorXd p_values(z_scores.size());
    boost::math::chi_squared chi_squared_dist(1);
    
    for (int i = 0; i < z_scores.size(); i++) {
        double chi_squared_stat = std::pow(z_scores(i), 2);
        p_values(i) = boost::math::cdf(boost::math::complement(chi_squared_dist, chi_squared_stat));
    }
    
    return p_values;
}

// Function to fit a Linear Mixed Model (LMM) for binary GWAS
std::tuple<std::string, std::string, std::string> LMM_binary(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, int>& binary_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covariate) {

    // Convert phenotype data into an Eigen vector
    Eigen::VectorXd y(binary_phenotype.size());
    Eigen::MatrixXd covariate_matrix(binary_phenotype.size(), covariate.begin()->second.size());
    int index = 0;
    std::vector<std::string> sample_ids;

    for (const auto& pair : binary_phenotype) {
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

    // Estimate beta using logistic regression approximation
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(design_matrix.cols());
    Eigen::VectorXd p = (1.0 / (1.0 + (-design_matrix * beta).array().exp())).matrix();
    Eigen::VectorXd W = p.array() * (1 - p.array());
    Eigen::MatrixXd XtWX = design_matrix.transpose() * W.asDiagonal() * design_matrix;
    Eigen::VectorXd XtWz = design_matrix.transpose() * (W.asDiagonal() * (design_matrix * beta + (y - p).array()).matrix());
    beta = XtWX.ldlt().solve(XtWz);

    // Compute p-values using Wald test
    Eigen::VectorXd p_values = computePValues(beta, XtWX);

    // Compute residuals
    Eigen::VectorXd residuals = y - (1.0 / (1.0 + (-design_matrix * beta).array().exp())).matrix();
    double residual_variance = (residuals.transpose() * residuals).value() / residuals.size();

    // Convert results to strings
    std::ostringstream beta_stream, p_values_stream, residual_var_stream;
    beta_stream << beta.transpose();
    p_values_stream << p_values.transpose();
    residual_var_stream << residual_variance;

    return std::make_tuple(beta_stream.str(), p_values_stream.str(), residual_var_stream.str());
}

// ------------------------ Chi2 test ------------------------

// Check if the observed matrix is valid (no zero rows/columns)
bool check_observed(const std::vector<std::vector<int>>& observed, size_t rows, size_t cols) {
    std::vector<int> col_sums(cols, 0);

    if (observed.size() == 0) return false;
    
    for (size_t i = 0; i < rows; ++i) {
        int row_sum = 0;
        for (size_t j = 0; j < cols; ++j) {
            row_sum += observed[i][j];
            col_sums[j] += observed[i][j];
        }
        if (row_sum <= 0) return false; // Check row sum
    }

    for (size_t j = 0; j < cols; ++j) {
        if (col_sums[j] <= 0) return false; // Check column sums
    }

    return true;
}

// Function to calculate the Chi-square test statistic
std::string chi2Test(const std::vector<std::vector<int>>& observed) {
    size_t rows = observed.size();
    size_t cols = observed[0].size();

    // Validate the observed matrix
    if (!check_observed(observed, rows, cols)) {
        return "NA";
    }

    // Compute row and column sums
    std::vector<double> row_sums(rows, 0.0);
    std::vector<double> col_sums(cols, 0.0);
    double total_sum = 0.0;

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            row_sums[i] += observed[i][j];
            col_sums[j] += observed[i][j];
            total_sum += observed[i][j];
        }
    }

    // Compute expected frequencies
    std::vector<std::vector<double>> expected(rows, std::vector<double>(cols, 0.0));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            expected[i][j] = (row_sums[i] * col_sums[j]) / total_sum;
        }
    }

    // Compute chi-squared statistic
    double chi_squared_stat = 0.0;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (expected[i][j] > 0) { // Avoid division by zero
                double diff = observed[i][j] - expected[i][j];
                chi_squared_stat += diff * diff / expected[i][j];
            }
        }
    }

    size_t degrees_of_freedom = (rows - 1) * (cols - 1);

    // Compute p-value using Boost's chi-squared distribution
    boost::math::chi_squared chi_squared_dist(degrees_of_freedom);
    double p_value = boost::math::cdf(boost::math::complement(chi_squared_dist, chi_squared_stat));

    return set_precision(p_value);
}

// ------------------------ Fisher exact test ------------------------

// Function to initialize the log factorials array
void initLogFacs(long double* logFacs, int n) {
    logFacs[0] = 0; 
    for (int i = 1; i < n+1; ++i) {
        logFacs[i] = logFacs[i - 1] + log((double)i);
    }
}

long double logHypergeometricProb(long double* logFacs , int a, int b, int c, int d) {
    return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d]
    - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

std::string fastFishersExactTest(const std::vector<std::vector<int>>& table) {
    // Ensure the table is 2x2
    if (table.size() != 2 || table[0].size() != 2 || table[1].size() != 2) {
        return "NA";
    }

    // Extract values from the table
    int a = table[0][0];
    int b = table[0][1];
    int c = table[1][0];
    int d = table[1][1];

    // Total sum of the table
    int n = a + b + c + d;
    long double* logFacs = new long double[n+1]; // *** dynamically allocate memory logFacs[0..n] ***
    initLogFacs(logFacs , n);

    long double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d);
    long double pFraction = 0;
    for(int x=0; x <= n; ++x) { // among all possible x
        int abx = a + b - x;
        int acx = a + c - x;
        int dax = d - a + x;
        if (abx >= 0 && acx >= 0 && dax >=0) { 
            long double l = logHypergeometricProb(logFacs, x, abx, acx, dax);
            if (l <= logpCutoff) {pFraction += exp(l - logpCutoff);}
        }
    }

    long double logpValue = exp(logpCutoff + log(pFraction));
    delete [] logFacs;
    return set_precision(logpValue);
}

// ------------------------ Binary table & stats ------------------------


std::string format_group_paths(const std::vector<std::vector<int>>& matrix) {

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

std::vector<std::string> binary_stat_test(const std::vector<std::vector<int>>& df) {

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
        int col_min = INT_MAX;
        for (const auto& row : df) {
            col_min = std::min(col_min, row[col]);
        }
        inter_group += col_min;
    }
    
    int average = static_cast<double>(allele_number) / numb_colum; // get 200 instead of 200.00000

    // Compute  Fisher's exact & Chi-squared test p-value
    std::string chi2_p_value = chi2Test(df);
    std::string fastfisher_p_value = fastFishersExactTest(df);
    std::string group_paths = format_group_paths(df); // Placeholder for future implementation

    return {fastfisher_p_value, chi2_p_value, std::to_string(allele_number), std::to_string(min_row_index), std::to_string(numb_colum), std::to_string(inter_group), std::to_string(average), group_paths};
}

std::vector<std::vector<int>> create_binary_table(
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, Matrix& matrix) 
{
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();
    size_t length_column_headers = list_path_snarl.size();

    // Initialize g0 and g1 with zeros, corresponding to the length of column_headers
    std::vector<int> g0(length_column_headers, 0);
    std::vector<int> g1(length_column_headers, 0);

    // Iterate over each path_snarl in column_headers
    for (size_t idx_g = 0; idx_g < list_path_snarl.size(); ++idx_g) {
        const std::string& path_snarl = list_path_snarl[idx_g];
        const size_t number_sample = list_samples.size();
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
