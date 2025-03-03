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

#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/distributions/students_t.hpp>

#include "matrix.hpp"
#include "snarl_parser.hpp"

// Linear regression function that returns a tuple of p_value, standard error (se), and beta
std::tuple<double, double, std::string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype);

std::unordered_map<std::string, std::vector<int>> create_quantitative_table(
    const std::vector<std::string>& column_headers,
    const std::vector<std::string>& list_samples, 
    Matrix& matrix);

#endif
