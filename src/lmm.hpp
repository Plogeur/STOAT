#ifndef LMM_H
#define LMM_H

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

#include "utils.hpp"
#include "arg_parser.hpp"

double normal_cdf(double x);

std::tuple<string, string, string, string> lmm_quantitative(
    const std::unordered_map<std::string, std::vector<size_t>>& df,                  
    const std::unordered_map<std::string, double>& phenotype_table,      
    const KinshipMatrix& kinship,                                              
    const std::unordered_map<std::string, std::vector<double>>& covariates);

#endif 