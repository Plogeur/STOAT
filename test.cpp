#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

double chi2Fast(const std::vector<std::vector<size_t>>& observed) {
    
    size_t cols = observed[0].size();
    std::vector<double> col_sums(cols, 0.0);
    double row_sum0 = 0.0, row_sum1 = 0.0, total = 0.0;

    // Precompute row sums and column sums
    for (size_t j = 0; j < cols; ++j) {
        double a = observed[0][j];
        double b = observed[1][j];
        row_sum0 += a;
        row_sum1 += b;
        col_sums[j] = a + b;
        total += a + b;
    }

    if (total == 0.0) return 0.0;

    // Compute chi-squared statistic
    double chi2 = 0.0;
    for (size_t j = 0; j < cols; ++j) {
        double expected0 = (row_sum0 * col_sums[j]) / total;
        double expected1 = (row_sum1 * col_sums[j]) / total;

        double diff0 = observed[0][j] - expected0;
        double diff1 = observed[1][j] - expected1;

        if (expected0 > 0) chi2 += (diff0 * diff0) / expected0;
        if (expected1 > 0) chi2 += (diff1 * diff1) / expected1;
    }

    size_t df = (observed.size() - 1) * (observed[0].size() - 1);

    boost::math::chi_squared dist(df);
    return boost::math::cdf(boost::math::complement(dist, chi2));
}

// Function to calculate the Chi-square test statistic
double chi2Test(const std::vector<std::vector<size_t>>& observed) {
    size_t rows = observed.size();
    size_t cols = observed[0].size();

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
    boost::math::chi_squared chi_squared_dist(degrees_of_freedom);
    return boost::math::cdf(boost::math::complement(chi_squared_dist, chi_squared_stat));
}

int main() {

    std::vector<std::vector<size_t>> observed = {
        {0, 0, 0, 0},
        {0, 40, 62, 10}
    };

    double pval_fast = 0.0;
    double pval_general = 0.0;

    // Time chi22_eval loop
    auto start_fast = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000; ++i) {
        pval_fast = chi2Fast(observed); // correct the chi22_eval to handle the observed matrix
    }
    auto end_fast = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_fast = end_fast - start_fast;

    // Time general chi2Test loop
    auto start_general = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100000; ++i) {
        pval_general = chi2Test(observed);
    }

    auto end_general = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_general = end_general - start_general;

    // Output results
    std::cout << "After 10,000 iterations:\n";
    std::cout << "Fast chi2Fast p-value: " << pval_fast << "\n";
    std::cout << "Time taken: " << duration_fast.count() << " seconds\n\n";

    std::cout << "General chi2Test p-value: " << pval_general << "\n";
    std::cout << "Time taken: " << duration_general.count() << " seconds\n";
    return 0;
}
