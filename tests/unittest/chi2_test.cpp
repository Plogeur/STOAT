#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <numeric>
#include <boost/math/distributions/chi_squared.hpp>

std::string set_precision(double value) {
    std::ostringstream oss;
    if (value < 0.0001) {
        oss << std::scientific << std::setprecision(4) << value; // Scientific notation with 4 decimals
    } else {
        oss << std::fixed << std::setprecision(4) << value; // Fixed-point notation with 4 decimals
    }
    return oss.str();
}

std::string chi2_2x2(const std::vector<int>& g0, const std::vector<int>& g1) {
    
    bool yates_correction = true;
    if (g0.size() != 2 || g1.size() != 2) {
        throw std::invalid_argument("Both input vectors must be of size 2 for a 2x2 table.");
    }

    int a = g0[0];
    int b = g0[1];
    int c = g1[0];
    int d = g1[1];

    int row1 = a + b;
    int row2 = c + d;
    int col1 = a + c;
    int col2 = b + d;
    int total = row1 + row2;

    if (row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0) {
        throw std::invalid_argument("Zero row or column detected in 2x2 table.");
    }

    double numerator = static_cast<double>(a * d - b * c);
    if (yates_correction) {
        numerator = std::abs(numerator) - 0.5 * total;
        numerator = std::max(0.0, numerator); // Prevent negative square root
    }

    numerator *= numerator;
    double denominator = static_cast<double>(row1 * row2 * col1 * col2) / total;

    double chi2_stat = numerator / denominator;

    // Get p-value using chi-squared distribution with 1 degree of freedom
    boost::math::chi_squared dist(1);
    double p_value = 1.0 - boost::math::cdf(dist, chi2_stat);

    return set_precision(p_value);
}

std::string chi2_test(const std::vector<int>& g0, const std::vector<int>& g1) {

    int cols = g0.size();
    std::vector<int> col_totals(cols);
    int total = 0;
    int row_total_0 = 0;
    int row_total_1 = 0;

    for (int i = 0; i < cols; ++i) {
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
    for (int i = 0; i < cols; ++i) {
        double expected_0 = static_cast<double>(row_total_0) * col_totals[i] / total;
        double expected_1 = static_cast<double>(row_total_1) * col_totals[i] / total;

        chi2 += (g0[i] - expected_0) * (g0[i] - expected_0) / expected_0;
        chi2 += (g1[i] - expected_1) * (g1[i] - expected_1) / expected_1;
    }

    int df = cols - 1;
    boost::math::chi_squared dist(df);
    return set_precision(1.0 - boost::math::cdf(dist, chi2));
}

void print_test(const std::vector<int>& g0, const std::vector<int>& g1, const std::string& name) {
    std::cout << "\n" << name << "\n";
    std::cout << "Group 0: ";
    for (auto v : g0) std::cout << std::setw(4) << v;
    std::cout << "\nGroup 1: ";
    for (auto v : g1) std::cout << std::setw(4) << v;
    std::cout << "\n";

    try {
        std::string p1 = chi2_test(g0, g1);
        std::cout << "Chi-squared test p-value: " << p1 << "\n";
        
        std::string p2 = chi2_2x2(g0, g1);
        std::cout << "Chi-squared 2X2 test p-value: " << p2 << "\n";

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }
}

int main() {
    std::vector<std::pair<std::vector<int>, std::vector<int>>> test_cases = {
        {{10, 20}, {20, 10}},                  // Balanced
        {{30, 5}, {2, 25}},                    // Strong effect
        {{10, 15, 5}, {20, 10, 10}},           // 2x3
        {{5, 10, 15, 20}, {20, 15, 10, 5}},    // 2x4
        {{10, 10, 10, 10, 10}, {10, 10, 10, 10, 10}}, // Uniform
        {{0, 0}, {0, 0}},                      // All Zeros
        {{0, 0, 0}, {10, 20, 30}},             // Full Zero Row
        {{0, 10, 5}, {0, 20, 15}},             // Full Zero Column
        {{0, 0}, {0, 1}},                      // One non-zero cell
        {{1, 0}, {0, 1}},                      // One non-zero cell
    };

    std::vector<std::string> names = {
        "Example 1: Balanced",
        "Example 2: Strong effect",
        "Example 3: 2x3 Table",
        "Example 4: 2x4 Table",
        "Example 5: All Zeros",
        "Example 6: Full Zero Row",
        "Example 7: Full Zero Column",
        "Example 8: Uniform Table",
        "Example 9: One Cell Non-Zero"
    };

    for (size_t i = 0; i < test_cases.size(); ++i) {
        print_test(test_cases[i].first, test_cases[i].second, names[i]);
    }

    return 0;
}
