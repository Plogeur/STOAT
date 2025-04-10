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
#include <Eigen/Dense>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

#include "matrix.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

// ------------------------ LMM BINARY ------------------------

std::vector<std::string> LMM_binary(
    const std::vector<std::vector<size_t>>& df,
    const std::unordered_map<std::string, std::vector<double>>& covariate);

// ------------------------ Chi2 test ------------------------

// Function to perform the Chi-square test
std::string chi2Test(const std::vector<std::vector<size_t>>& observed);

// ------------------------ Fisher exact test ------------------------

// Function to perform Fisher's exact test
std::string fastFishersExactTest(const std::vector<std::vector<size_t>>& table);

// ------------------------ Binary table ------------------------

std::string format_group_paths(const std::vector<std::vector<size_t>>& df);

std::vector<std::string> binary_stat_test(const std::vector<std::vector<size_t>>& df);

std::vector<std::vector<size_t>> create_binary_table(
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, const Matrix& matrix);

#endif
