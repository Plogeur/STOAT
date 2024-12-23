#ifndef QUANTITATIVE_ANALYSIS_HPP
#define QUANTITATIVE_ANALYSIS_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include "matrix.hpp"
#include "snarl_parser.hpp"

// Function to calculate the mean of a vector
double mean(const std::vector<double>& v);

// Function to calculate the variance of a vector
double variance(const std::vector<double>& v);

// Function to calculate the covariance between two vectors
double covariance(const std::vector<double>& x, const std::vector<double>& y);

// Linear regression function that returns a tuple of p_value, standard error (se), and beta
std::tuple<double, double, double> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, float>& quantitative_phenotype);

std::unordered_map<std::string, std::vector<int>> create_quantitative_table(
    const std::vector<std::string>& column_headers,
    const std::vector<std::string>& list_samples, 
    Matrix& matrix);

#endif
