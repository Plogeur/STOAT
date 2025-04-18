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

// ------------------------ Logistic regression ------------------------

double sigmoid(double z);
double dot(const std::vector<double>& a, const std::vector<double>& b);
double normal_p_value(double z);

std::tuple<string, string, string, string> logistic_regression(
    const std::unordered_map<std::string, std::vector<size_t>>& variant_data,
    const std::unordered_map<std::string, bool>& phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covariates,
    double tol, size_t max_iter);

// ------------------------ Chi2 test ------------------------

// Function to perform the Chi-square test
std::string chi2Test(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// ------------------------ Fisher exact test ------------------------

// Function to perform Fisher's exact test
std::string fastFishersExactTest(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// ------------------------ Binary table ------------------------

std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

void binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1, 
    string& fastfisher_p_value, string& chi2_p_value, string& group_paths,
    string& allele_number_str, string& min_row_index_str, string& numb_colum_str, 
    string& inter_group_str, string& average_str);

bool create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl, 
    const std::vector<std::string>& list_samples, 
    const Matrix& matrix, const double& maf);

#endif
