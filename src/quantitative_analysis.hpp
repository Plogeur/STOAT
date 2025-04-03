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
std::tuple<string, string, string, string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype, 
    const size_t& total_snarl);

std::pair<std::unordered_map<std::string, std::vector<int>>, size_t> create_quantitative_table(
    const std::vector<std::string>& list_samples, 
    const std::vector<std::string>& column_headers,
    Matrix& matrix);

std::tuple<std::string, std::string, std::string> LMM_quantitative(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype,
    const std::unordered_map<std::string, std::vector<double>>& covariate);

#endif
