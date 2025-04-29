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
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "matrix.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

// ------------------------ Logistic regression ------------------------

double sigmoid(double z);
double compute_r2(const Eigen::VectorXd& y, const Eigen::VectorXd& p_null, const Eigen::VectorXd& p_full);

void logistic_regression(
    const std::vector<std::vector<size_t>>& variant_data,
    const std::vector<bool>& phenotype,
    const std::vector<string>& list_samples,
    const std::unordered_map<std::string, std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str);

// ------------------------ Chi2 test ------------------------

// Function to perform the Chi-square test on row size > 2 
std::string chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// Function to perform the Chi-square test on row size == 2 
std::string chi2_2x2(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// ------------------------ Fisher exact test ------------------------

// Function to perform Fisher's exact test
std::string fastFishersExactTest(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// ------------------------ Binary table ------------------------

std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

void binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1, 
    string& fastfisher_p_value, string& chi2_p_value, string& group_paths,
    string& allele_number_str, string& min_row_index_str, string& numb_colum_str, 
    string& inter_group_str, string& average_str);

size_t create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const vector<bool>& binary_phenotype, 
    const std::vector<std::string>& list_path_snarl, 
    const size_t& number_paths,
    const size_t& number_samples,
    const Matrix& matrix);

bool check_MAF_threshold(
    const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    const size_t& totalSum, const size_t& length_column_headers, 
    const double& maf);

#endif
