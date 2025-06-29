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
#include <Eigen/Core>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "matrix.hpp"
#include "snarl_analyser.hpp"
#include "utils.hpp"

// ------------------------ Logistic regression ------------------------
double normal_cdf(double z);
inline double sigmoid(double x);
inline double clamp(double x, double lo, double hi);
double calculate_log_likelihood(const Eigen::VectorXd& y, const Eigen::VectorXd& p);

// Standard normal cumulative distribution function
double normal_cdf(double z);

// Sigmoid function
inline double sigmoid(double x);

// Clamp helper
inline double clamp(double x, double lo, double hi);

void logistic_regression(
    const std::vector<std::vector<double>>& variants_data,
    const std::vector<bool>& phenotype,
    std::string& p_value_out,
    std::string& beta_out,
    std::string& se_out,
    std::string& r2_out);

void glm_logistic_covar(
    const std::vector<std::vector<double>>& variant_data,
    const std::vector<bool>& phenotype,
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str);
    
// ------------------------ Chi2 test ------------------------

// Function to perform the Chi-square test on row size > 2 
std::string chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// Function to perform the Chi-square test on row size == 2 
std::string chi2_2x2(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// ------------------------ Fisher exact test ------------------------

// Function to perform Fisher's exact test
std::string fastFishersExactTest(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// ------------------------ Binary table ------------------------

// Write a string of: g0[0]:g1[1],g0[1]:g1[1],g0[2]:g1[2]...
std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// Given two vectors of genotypes representing two groups, fill in the p-values, etc by running the relevant tests
void binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1, 
    string& fastfisher_p_value, string& chi2_p_value, string& group_paths,
    string& allele_number_str, string& min_row_index_str, string& numb_colum_str, 
    string& inter_group_str, string& average_str);

// Given two vectors of genotypes representing two groups (with length number_paths), fill them in with counts of the number of times each path is seen  
// g0 and g1 can be used in binary_stat_test()
size_t create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const vector<bool>& binary_phenotype, 
    const std::vector<Path_traversal_t>& list_path_snarl, 
    const size_t& number_paths,
    const size_t& number_samples,
    const EdgeBySampleMatrix& matrix);

// Does at least one column meet the MAF threshold?
bool check_MAF_threshold(
    const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    const size_t& totalSum, const size_t& length_column_headers, 
    const double& maf);

#endif
