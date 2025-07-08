#include "binary_table.hpp"

// ------------------------ Binary table & stats ------------------------
namespace stoat_vcf {

std::tuple<std::string, 
std::string, std::string, std::string, std::string, 
std::string> binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1) {

    // Compute derived statistics
    int allele_number = 0;
    int inter_group = 0;
    int numb_colum = g0.size();
    int min_row_index = INT_MAX;

    std::string
    group_paths, allele_number_str, min_row_index_str, 
    numb_colum_str, inter_group_str, average_str;

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

    group_paths = stoat_vcf::format_group_paths(g0, g1);
    allele_number_str = std::to_string(allele_number);
    min_row_index_str = std::to_string(min_row_index);
    numb_colum_str = std::to_string(numb_colum);
    inter_group_str = std::to_string(inter_group);
    average_str = std::to_string(average);

    return std::make_tuple(
        group_paths, allele_number_str, min_row_index_str, 
        numb_colum_str, inter_group_str, average_str);
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

} // namespace stoat_vcf
