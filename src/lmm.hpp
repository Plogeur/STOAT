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

void lmm_quantitative(
    const std::vector<std::vector<size_t>>& df,                  
    const vector<double>& phenotype_table,      
    const KinshipMatrix& kinship,                                              
    const std::unordered_map<std::string, std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str);

void lmm_binary(
    const std::vector<std::vector<size_t>>& df,                  
    const std::vector<bool>& phenotype_binary,     
    const KinshipMatrix& kinship,                                              
    const std::unordered_map<std::string, std::vector<double>>& covariates,
    std::string& p_value_str, std::string& beta_str, std::string& se_str, std::string& r2_str);

#endif 