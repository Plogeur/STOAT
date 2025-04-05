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

using namespace std;

std::string set_precision(double value);
double mean_pvalue_from_strings(const std::string& p1, const std::string& p2);
bool is_na(const std::string& s);
double string_to_pvalue(const std::string& p1);

#endif
