#include "binary_test.hpp"
#include "snarl_analyser.hpp"
#include "utils.hpp"

// ------------------------ Binary table & stats ------------------------
namespace stoat_vcf {

void binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    std::string& fastfisher_p_value, std::string& chi2_p_value, std::string& group_paths,
    std::string& allele_number_str, std::string& min_row_index_str, std::string& numb_colum_str, 
    std::string& inter_group_str, std::string& average_str) {

    // Compute derived statistics
    int allele_number = 0;
    int inter_group = 0;
    int numb_colum = g0.size();
    int min_row_index = INT_MAX;

    for (size_t i = 0; i < g0.size(); ++i) {
        int row_sum = static_cast<int>(g0[i] + g1[i]);
        allele_number += row_sum;
        min_row_index = std::min(min_row_index, row_sum);
    }

    for (int col=0; col < numb_colum; ++col) {
        size_t col_min = INT_MAX;
        col_min = std::min(col_min, g0[col]);
        col_min = std::min(col_min, g1[col]);
        inter_group += col_min;
    }
    
    int average = static_cast<double>(allele_number) / numb_colum; // get 200 instead of 200.00000

    // Compute  Fisher's exact & Chi-squared test p-value
    if (g0.size() == 2) {
        size_t a = g0[0];
        size_t b = g0[1];
        size_t c = g1[0];
        size_t d = g1[1];
        chi2_p_value = chi2_2x2(a, b, c, d);
        fastfisher_p_value = fastFishersExactTest(a, b, c, d);
    } else {
        chi2_p_value = chi2_2xN(g0, g1);
    }
    group_paths = stoat_vcf::format_group_paths(g0, g1);
    allele_number_str = std::to_string(allele_number);
    min_row_index_str = std::to_string(min_row_index);
    numb_colum_str = std::to_string(numb_colum);
    inter_group_str = std::to_string(inter_group);
    average_str = std::to_string(average);
}

size_t create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const std::vector<bool>& binary_phenotype, 
    const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl, 
    const size_t& number_paths,
    const size_t& number_samples,
    const stoat_vcf::EdgeBySampleMatrix& matrix) {

    size_t total_sum = 0;
    for (size_t idx_g = 0; idx_g < number_paths; ++idx_g) {
        const stoat_vcf::Path_traversal_t& path_snarl = list_path_snarl[idx_g];
        std::vector<stoat_vcf::Edge_t> list_edge_path = stoat_vcf::decompose_path_to_edges(path_snarl);
        std::vector<size_t> idx_srr_save = stoat_vcf::identify_path(list_edge_path, matrix, number_samples * 2);

        for (size_t idx : idx_srr_save) {
            bool group = binary_phenotype[idx / 2];
            if (group) {
                g1[idx_g] += 1;
            } else {
                g0[idx_g] += 1;
            }
            total_sum++;
        }
    }
    return total_sum;
}

bool check_MAF_threshold(
    const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    const size_t& totalSum, const size_t& length_column_headers, 
    const double& maf) {

    // Check MAF threshold
    for (size_t i = 0; i < length_column_headers; ++i) {
        int columnSum = g0[i] + g1[i];
        if (static_cast<double>(columnSum) / totalSum >= maf) {
            return true; // MAF threshold met
        }
    }
    return false; // No column met MAF threshold
}

} // namespace stoat_vcf
