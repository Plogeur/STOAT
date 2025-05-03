#ifndef QUANTITATIVE_ANALYSIS_HPP
#define QUANTITATIVE_ANALYSIS_HPP

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
#include <boost/math/distributions/chi_squared.hpp>

#include "matrix.hpp"
#include "snarl_parser.hpp"
#include "utils.hpp"

// Linear regression function that returns a tuple of p_value, standard error (se), and beta
void linear_regression(
    const std::vector<std::vector<size_t>>& df,
    const std::vector<double>& quantitative_phenotype,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str);

void glm_quantitative(
    const std::vector<std::vector<size_t>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str);

std::pair<std::vector<std::vector<size_t>>, size_t> create_quantitative_table(
    const size_t& length_sample,
    const std::vector<std::string>& column_headers,
    Matrix& matrix);

#endif
