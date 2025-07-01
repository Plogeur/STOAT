#ifndef stats_test_HPP
#define stats_test_HPP

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

#include "arg_parser.hpp"
#include "matrix.hpp"
#include "snarl_analyser.hpp"
#include "utils.hpp"

using namespace std;

// ------------------------ Linear regression ------------------------

// Linear regression function OLS with intercept + covariate if not empty

void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar,
    std::string& p_value_str, 
    std::string& beta_str, 
    std::string& se_str, 
    std::string& r2_str);

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
    const std::vector<std::vector<double>>& variant_data,
    const std::vector<bool>& phenotype,
    const std::vector<std::vector<double>>& covariates,
    std::string& p_value_str, 
    std::string& beta_str, 
    std::string& se_str, 
    std::string& r2_str);
    
// ------------------------ Chi2 test ------------------------

// Function to perform the Chi-square test on row size > 2 
std::string chi2_2xN(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// Function to perform the Chi-square test on row size == 2 
std::string chi2_2x2(const size_t& m11, const size_t& m12,
    const size_t& m21, const size_t& m22);

// ------------------------ Fisher exact test ------------------------

// Function to perform Fisher's exact test
std::string fastFishersExactTest(size_t m11, size_t m12,
    size_t m21, size_t m22);

// ------------------------------ LMM ------------------------------

// void lmm_quantitative(
//     const std::vector<std::vector<double>>& df,
//     const std::vector<double>& phenotype_table,
//     const stoat_vcf::KinshipMatrix& kinship,
//     const std::vector<std::vector<double>>& covariates,
//     std::string& p_value_str, 
//     std::string& beta_str, 
//     std::string& se_str, 
//     std::string& r2_str);

// void lmm_binary(
//     const std::vector<std::vector<double>>& df,
//     const std::vector<bool>& phenotype_binary,
//     const stoat_vcf::KinshipMatrix& kinship,
//     const std::vector<std::vector<double>>& covariates,
//     std::string& p_value_str, 
//     std::string& beta_str, 
//     std::string& se_str, 
//     std::string& r2_str);

#endif 
