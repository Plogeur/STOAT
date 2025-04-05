#ifndef ADJUSTED_PVALUE_HPP
#define ADJUSTED_PVALUE_HPP

#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>

void adjust_pvalues_BH(
    std::vector<std::tuple<double, double, size_t>>& data);

void update_gwas_file_with_adjusted_pvalues(
    const std::string& input_file,
    const std::string& output_file_significant,
    const int column_index,
    const std::vector<std::tuple<double, double, size_t>>& adjusted_pvalues);


#endif // ADJUSTED_PVALUE_HPP