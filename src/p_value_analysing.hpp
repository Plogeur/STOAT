#ifndef P_VALUE_ANALYSING_HPP
#define P_VALUE_ANALYSING_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>

struct Variant_Binary {
    std::string chr, pos, snarl, type;
    double p_fisher, p_chi2, p_adjusted = 0.0;
    int allele_num, min_row_index, num_column;
    std::string inter_group, average, group_paths;
};

bool compare_PValue_Binary(const Variant_Binary &a, const Variant_Binary &b);
void benjamini_Hochberg_Correction_Binary(std::vector<Variant_Binary> &variants);
void process_Binary(const std::string &input_file, const std::string &top_variant_file);

struct Variant_Quantitative {
    std::string chr, pos, snarls, type, allele_num;
    double p_value;
    double p_adjusted = 0.0;
};

bool compare_PValue_Quantitative(const Variant_Quantitative &a, const Variant_Quantitative &b);
void benjamini_Hochberg_Correction_Quantitative(std::vector<Variant_Quantitative> &variants);
void process_Quantitative(const std::string &input_file, const std::string &top_variant_file);

#endif
