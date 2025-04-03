#include "p_value_analysing.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <limits>

// Compare the mean of the Fisher and Chi2 p-values for sorting
bool compare_PValue_Binary(const Variant_Binary &a, const Variant_Binary &b) {
    // Calculate the mean p-value while handling NA
    double mean_pvalue_a = (std::isnan(a.p_fisher) ? a.p_chi2 : (std::isnan(a.p_chi2) ? a.p_fisher : (a.p_fisher + a.p_chi2) / 2));
    double mean_pvalue_b = (std::isnan(b.p_fisher) ? b.p_chi2 : (std::isnan(b.p_chi2) ? b.p_fisher : (b.p_fisher + b.p_chi2) / 2));

    return mean_pvalue_a < mean_pvalue_b; // Sorting based on mean p-value for BH correction
}

// Apply Benjamini-Hochberg correction based on the mean p-value
void benjamini_Hochberg_Correction_Binary(std::vector<Variant_Binary> &variants) {
    size_t m = variants.size();
    std::vector<size_t> indices(m);
    for (size_t i = 0; i < m; ++i) indices[i] = i;
    
    // Sort by mean p-value (based on Fisher and Chi2 combined), considering NA values
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        double mean_pvalue_i = (std::isnan(variants[i].p_fisher) ? variants[i].p_chi2 : (std::isnan(variants[i].p_chi2) ? variants[i].p_fisher : (variants[i].p_fisher + variants[i].p_chi2) / 2));
        double mean_pvalue_j = (std::isnan(variants[j].p_fisher) ? variants[j].p_chi2 : (std::isnan(variants[j].p_chi2) ? variants[j].p_fisher : (variants[j].p_fisher + variants[j].p_chi2) / 2));
        return mean_pvalue_i < mean_pvalue_j;
    });

    // Apply BH correction
    for (size_t i = 0; i < m; ++i) {
        size_t idx = indices[i];
        double rank = i + 1;
        double mean_pvalue = (std::isnan(variants[idx].p_fisher) ? variants[idx].p_chi2 : (std::isnan(variants[idx].p_chi2) ? variants[idx].p_fisher : (variants[idx].p_fisher + variants[idx].p_chi2) / 2));
        variants[idx].p_adjusted = std::min(1.0, mean_pvalue * m / rank);
    }
}

// Process the binary data file and apply the corrections
void process_Binary(const std::string &input_file, const std::string &top_variant_file) {
    std::ifstream infile(input_file);
    std::ofstream topfile(top_variant_file);
    
    if (!infile || !topfile) {
        std::cerr << "Error opening files." << std::endl;
        return;
    }
    
    std::vector<Variant_Binary> variants;
    std::string line;
    getline(infile, line); // Read and write header
    topfile << line << "\n";
    
    while (getline(infile, line)) {
        std::istringstream iss(line);
        Variant_Binary v;
        
        // Read the values and handle NA (use NaN for NA values)
        iss >> v.chr >> v.pos >> v.snarl >> v.type 
            >> v.p_fisher >> v.p_chi2 >> v.allele_num 
            >> v.min_row_index >> v.num_column
            >> v.inter_group >> v.average >> v.group_paths;
        
        // Convert missing p-values ("NA") into NaN (Not a Number) for calculation
        if (v.p_fisher == -1) v.p_fisher = std::numeric_limits<double>::quiet_NaN();
        if (v.p_chi2 == -1) v.p_chi2 = std::numeric_limits<double>::quiet_NaN();
        
        variants.push_back(v);
    }
    
    // Apply Benjamini-Hochberg correction
    benjamini_Hochberg_Correction_Binary(variants);

    // Open the input file again in output mode to overwrite the file
    infile.close();
    std::ofstream outfile(input_file);

    // Write the updated file with adjusted p-values
    outfile << "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tP_ADJUSTED\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n";
    for (const auto &v : variants) {
        outfile << v.chr << "\t" << v.pos << "\t" << v.snarl << "\t" << v.type << "\t"
                << v.p_fisher << "\t" << v.p_chi2 << "\t" << v.p_adjusted << "\t"
                << v.allele_num << "\t" << v.min_row_index << "\t"
                << v.num_column << "\t" << v.inter_group << "\t"
                << v.average << "\t" << v.group_paths << "\n";
        
        // Write top variants with adjusted p-value < 10^-5
        if (v.p_adjusted < 1e-5) {
            topfile << v.chr << "\t" << v.pos << "\t" << v.snarl << "\t" << v.type << "\t"
                    << v.p_fisher << "\t" << v.p_chi2 << "\t" << v.p_adjusted << "\t"
                    << v.allele_num << "\t" << v.min_row_index << "\t"
                    << v.num_column << "\t" << v.inter_group << "\t"
                    << v.average << "\t" << v.group_paths << "\n";
        }
    }
    
    infile.close();
    outfile.close();
    topfile.close();
}

bool compare_PValue_Quantitative(const Variant_Quantitative &a, const Variant_Quantitative &b) {
    return a.p_value < b.p_value; // Sorting based on p_fisher for BH correction
}

void benjamini_Hochberg_Correction_Quantitative(std::vector<Variant> &variants) {
    size_t m = variants.size();
    std::vector<size_t> indices(m);
    for (size_t i = 0; i < m; ++i) indices[i] = i;
    
    std::sort(indices.begin(), indices.end(), [&](size_t i, size_t j) {
        return variants[i].p_value < variants[j].p_value;
    });
    
    for (size_t i = 0; i < m; ++i) {
        size_t idx = indices[i];
        double rank = i + 1;
        variants[idx].p_adjusted = std::min(1.0, variants[idx].p_value * m / rank);
    }
}

void process_Quantitative(const std::string &input_file, const std::string &top_variant_file) {
    std::ifstream infile(input_file);
    std::ofstream topfile(top_variant_file);
    
    std::vector<Variant_Quantitative> variants;
    std::string line;
    getline(infile, line); // Read and write header
    topfile << line << "\n";
    
    while (getline(infile, line)) {
        std::istringstream iss(line);
        Variant_Quantitative v;
        iss >> v.chr >> v.pos >> v.snarl >> v.type >> v.beta >> v.se >> v.p_value >> v.p_adjusted >> v.allele_num;
        variants.push_back(v);
    }
    
    benjamini_Hochberg_Correction_Quantitative(variants);

    // Open the input file again in output mode to overwrite the file
    infile.close();
    std::ofstream outfile(input_file);
    
    // Write the updated file with adjusted p-values
    outfile << "CHR\tPOS\tSNARL\tTYPE\tBETA\tSE\tP\tP_ADJUSTED\tALLELE_NUM\n";
    for (const auto &v : variants) {
        outfile << v.chr << "\t" << v.pos << "\t" << v.snarl << "\t" << v.type << "\t"
                << v.beta << "\t" << v.se << "\t" << v.p_value << "\t" << v.p_adjusted << "\t"
                << v.allele_num << "\n";
        
        if (v.p_adjusted < 1e-5) {
            topfile << v.chr << "\t" << v.pos << "\t" << v.snarl << "\t" << v.type << "\t"
                    << v.beta << "\t" << v.se << "\t" << v.p_value << "\t" << v.p_adjusted << "\t"
                    << v.allele_num << "\n";
        }
    }
    
    infile.close();
    outfile.close();
    topfile.close();
}