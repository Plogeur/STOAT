#include "binary_analysis.hpp"
#include "snarl_parser.hpp"
#include <sstream>

// String utilities
std::string join(const std::vector<std::string>& elements, const std::string& delimiter) {
    if (elements.empty()) return "";
    
    std::ostringstream oss;
    auto it = elements.begin();
    oss << *it++;
    
    while (it != elements.end()) {
        oss << delimiter << *it++;
    }
    
    return oss.str();
}

// ------------------------ Chi2 test ------------------------

// Check if the observed matrix is valid (no zero rows/columns)
bool check_observed(const std::vector<std::vector<int>>& observed, size_t rows, size_t cols) {
    if (observed.empty() || observed[0].empty()) return false;
    
    // Check for zero rows
    for (size_t i = 0; i < rows; ++i) {
        int row_sum = 0;
        for (size_t j = 0; j < cols; ++j) {
            row_sum += observed[i][j];
        }
        if (row_sum == 0) return false;
    }
    
    // Check for zero columns
    for (size_t j = 0; j < cols; ++j) {
        int col_sum = 0;
        for (size_t i = 0; i < rows; ++i) {
            col_sum += observed[i][j];
        }
        if (col_sum == 0) return false;
    }
    
    return true;
}

// Function to calculate the Chi-square test statistic
std::string chi2Test(const std::vector<std::vector<int>>& observed) {
    if (observed.empty() || observed[0].empty()) return "NA";
    
    size_t rows = observed.size();
    size_t cols = observed[0].size();
    
    if (!check_observed(observed, rows, cols)) return "NA";
    
    // Calculate row and column totals
    std::vector<int> row_totals(rows, 0);
    std::vector<int> col_totals(cols, 0);
    int total = 0;
    
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            row_totals[i] += observed[i][j];
            col_totals[j] += observed[i][j];
            total += observed[i][j];
        }
    }
    
    // Calculate chi-square statistic
    double chi_square = 0.0;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            double expected = static_cast<double>(row_totals[i] * col_totals[j]) / total;
            double diff = observed[i][j] - expected;
            chi_square += (diff * diff) / expected;
        }
    }
    
    // Calculate degrees of freedom
    int df = (rows - 1) * (cols - 1);
    
    // Calculate p-value using boost
    try {
        boost::math::chi_squared dist(df);
        double p_value = 1.0 - boost::math::cdf(dist, chi_square);
        
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4) << p_value;
        return oss.str();
    } catch (...) {
        return "NA";
    }
}

// ------------------------ Fisher exact test ------------------------

// Function to initialize the log factorials array
void initLogFacs(long double* logFacs, int n) {
    logFacs[0] = 0.0;
    for (int i = 1; i <= n; ++i) {
        logFacs[i] = logFacs[i-1] + std::log(i);
    }
}

long double logHypergeometricProb(long double* logFacs, int a, int b, int c, int d) {
    if (a < 0 || b < 0 || c < 0 || d < 0) return 0.0;
    
    int r1 = a + b;
    int r2 = c + d;
    int c1 = a + c;
    int c2 = b + d;
    int n = r1 + r2;
    
    return logFacs[r1] + logFacs[r2] + logFacs[c1] + logFacs[c2] 
           - logFacs[n] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d];
}

long double fastFishersExactTest(const std::vector<std::vector<int>>& table) {
    if (table.size() != 2 || table[0].size() != 2) return -1.0;
    
    int a = table[0][0];
    int b = table[0][1];
    int c = table[1][0];
    int d = table[1][1];
    
    if (a < 0 || b < 0 || c < 0 || d < 0) return -1.0;
    if (a == 0 && b == 0 && c == 0 && d == 0) return -1.0;
    
    int n = a + b + c + d;
    long double logFacs[n + 1];
    initLogFacs(logFacs, n);
    
    long double observed = std::exp(logHypergeometricProb(logFacs, a, b, c, d));
    long double sum = 0.0;
    
    int min_a = std::max(0, a + b - d);
    int max_a = std::min(a + b, a + c);
    
    for (int i = min_a; i <= max_a; ++i) {
        int curr_b = a + b - i;
        int curr_c = a + c - i;
        int curr_d = n - i - curr_b - curr_c;
        
        long double prob = std::exp(logHypergeometricProb(logFacs, i, curr_b, curr_c, curr_d));
        if (prob <= observed) sum += prob;
    }
    
    return sum;
}

// ------------------------ Binary table & stats ------------------------


std::string format_group_paths(const std::vector<std::vector<int>>& df) {
    if (df.empty()) return "";
    
    std::vector<std::string> formatted;
    for (size_t j = 0; j < df[0].size(); ++j) {
        std::string col;
        for (size_t i = 0; i < df.size(); ++i) {
            if (i > 0) col += ":";
            col += std::to_string(df[i][j]);
        }
        formatted.push_back(col);
    }
    
    return join(formatted, ",");
}

std::vector<std::string> binary_stat_test(const std::vector<std::vector<int>>& df) {
    std::vector<std::string> results(8, "NA");
    
    if (df.empty() || df[0].empty()) {
        results[2] = "0";
        results[3] = std::to_string(std::numeric_limits<int>::max());
        results[4] = "0";
        results[5] = "0";
        results[6] = "0.000000";
        results[7] = "";
        return results;
    }
    
    std::vector<std::vector<int>> contingency_table(2, std::vector<int>(2, 0));
    contingency_table[0][0] = df[0][0];
    contingency_table[0][1] = df[0][1];
    contingency_table[1][0] = df[1][0];
    contingency_table[1][1] = df[1][1];
    
    results[0] = std::to_string(fastFishersExactTest(contingency_table));
    results[1] = chi2Test(contingency_table);
    results[2] = std::to_string(df[0].size());
    
    int min_row = std::numeric_limits<int>::max();
    for (const auto& row : df) {
        for (int val : row) {
            min_row = std::min(min_row, val);
        }
    }
    results[3] = std::to_string(min_row);
    
    results[4] = std::to_string(df[0].size());
    
    int inter_group = 0;
    for (size_t j = 0; j < df[0].size(); ++j) {
        bool has_both = false;
        for (size_t i = 0; i < df.size(); ++i) {
            if (df[i][j] > 0) {
                has_both = true;
                break;
            }
        }
        if (has_both) inter_group++;
    }
    results[5] = std::to_string(inter_group);
    
    double average = 0.0;
    int count = 0;
    for (const auto& row : df) {
        for (int val : row) {
            average += val;
            count++;
        }
    }
    average /= count;
    
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << average;
    results[6] = oss.str();
    
    results[7] = format_group_paths(df);
    
    return results;
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
