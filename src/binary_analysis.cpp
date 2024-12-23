#include "binary_analysis.hpp"
#include "snarl_parser.hpp"

// ------------------------ Chi2 test ------------------------

// Check if the observed matrix is valid (no zero rows/columns)
bool check_observed(const std::vector<std::vector<int>>& observed, size_t rows, size_t cols) {
    std::vector<int> col_sums(cols, 0);

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

    // Format p-value as a string
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(4) << p_value;
    return ss.str();
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

long double fastFishersExactTest(const std::vector<std::vector<int>>& table) {
    // Ensure the table is 2x2
    if (table.size() != 2 || table[0].size() != 2 || table[1].size() != 2) {
        throw std::invalid_argument("Input table must be 2x2.");
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
        if ( abx >= 0 && acx >= 0 && dax >=0 ) { 
            long double l = logHypergeometricProb(logFacs, x, abx, acx, dax);
            if (l <= logpCutoff) {pFraction += exp(l - logpCutoff);}
        }
    }

    long double logpValue = exp(logpCutoff + log(pFraction));
    delete [] logFacs;
    return logpValue;
}

// ------------------------ Binary table & stats ------------------------

std::vector<std::string> binary_stat_test(const std::vector<std::vector<int>>& df) {
    long double fastfisher_p_value = fastFishersExactTest(df);

    std::string chi2_p_value = chi2Test(df);
    int group_I_path_I = df[0][0];
    int group_I_path_II = df[0][1];
    int group_II_path_I = df[1][0];
    int group_II_path_II = df[1][1];

    // Pvalue Precision of 6 number after 0.
    std::stringstream ss;
    ss << std::scientific << std::setprecision(4) << fastfisher_p_value;
    std::string stringFastfisher_p_value = ss.str();

    std::vector<std::string> result = {
        stringFastfisher_p_value,
        chi2_p_value,
        std::to_string(group_I_path_I),
        std::to_string(group_I_path_II),
        std::to_string(group_II_path_I),
        std::to_string(group_II_path_II)
    };
    return result;
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
