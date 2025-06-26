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

// Given a vector of <p-value, 1.0, line index from the input file>, fill in the vector with the adjusted p-value
// and sort the vector by adjusted p-value
// TODO: This header isn't actually what got implemented
void adjust_pvalues_BH(std::vector<double>& vector_pvalues);

// Read a tsv from input_file, collect the p-values from the correct column (depending on phenotype_type), 
// and write the same file plus a BH-adjusted p-value to outupt_file_significant.
void add_BH_adjusted_column(
    const std::string& input_file, 
    const std::string& output_file_significant,
    const std::string& phenotype_type);

#endif // ADJUSTED_PVALUE_HPP
