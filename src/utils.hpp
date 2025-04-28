#ifndef UTILS_HPP
#define UTILS_HPP

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
#include <fstream>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;

std::string set_precision(const double& value);
std::string set_precision_chi2(const boost::multiprecision::cpp_dec_float_50& value);

double combine_pvalue_from_strings(const std::string& p1, const std::string& p2);
bool is_na(const std::string& s);
double string_to_pvalue(const std::string& p1);
void writeSignificantTableToTSV(
    const std::unordered_map<std::string, std::vector<size_t>>& table,
    const std::vector<std::string>& list_snarl,
    const std::string& filename);

bool isPValueSignificant(size_t numDigits, const std::string& pvalue_str);

#endif
