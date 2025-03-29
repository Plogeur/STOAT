#ifndef BINARY_ANALYSIS_HPP
#define BINARY_ANALYSIS_HPP

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <sstream>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

#include "matrix.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

// ------------------------ LMM ------------------------

void LMM_binary(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, int>& binary_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covariate);

// ------------------------ Chi2 test ------------------------

// Function to check if matrix is valid (no zero rows/columns)
bool check_observed(const std::vector<std::vector<int>>& observed, size_t rows, size_t cols);

// Function to perform the Chi-square test
std::string chi2Test(const std::vector<std::vector<int>>& observed);

// ------------------------ Fisher exact test ------------------------

long double logHypergeometricProb(long double* logFacs , int a, int b, int c, int d);
void initLogFacs(long double* logFacs, int n);

// Function to perform Fisher's exact test
std::string fastFishersExactTest(const std::vector<std::vector<int>>& table);

// ------------------------ Binary table ------------------------

std::string format_group_paths(const std::vector<std::vector<int>>& df);

std::vector<std::string> binary_stat_test(const std::vector<std::vector<int>>& df);

std::vector<std::vector<int>> create_binary_table(
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, Matrix& matrix);

#endif
