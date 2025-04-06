#ifndef POST_PROCESSING_HPP
#define POST_PROCESSING_HPP

#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>

void adjust_pvalues_BH(std::vector<double>& vector_pvalues);

void add_BH_adjusted_column(
    const std::string& input_file, 
    const std::string& output_file_significant,
    const std::string& phenotype_type);

#endif // ADJUSTED_PVALUE_HPP