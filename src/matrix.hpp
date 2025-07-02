#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

#include "snarl_data_t.hpp"

using namespace std;

namespace stoat_vcf {

// A class to store a 2d bit-matrix
// Rows represent edges and the index of each edge can be found from the row_header
// Columns represent samples/haplotypes
class EdgeBySampleMatrix {
public:
    EdgeBySampleMatrix(size_t rows, size_t cols);
    ~EdgeBySampleMatrix()=default;

    // Operator to get the value
    bool operator()(size_t row, size_t col) const;

    // Add this edge to the matrix
    void push_matrix(const Edge_t& EdgePath, std::unordered_map<stoat_vcf::Edge_t, size_t>& edge_dict, size_t indexColumn);

    // Set this value to true
    void set(size_t row, size_t col);

    // Get the matrix itself
    const std::vector<uint8_t>& get_matrix() const;

    // Get the row_header
    const std::unordered_map<stoat_vcf::Edge_t, size_t>& get_row_header() const;

    // Get the maximum index into the vector representing the matrix
    size_t getMaxElement() const;
    
    // Double the size of the matrix
    void expandMatrix();

    // Reset row_header
    void set_row_header(const std::unordered_map<stoat_vcf::Edge_t, size_t>& row_header);

    // Shrink to use the minimum amount of memory possible allowing current_rows
    void shrink(size_t current_rows);

    // Return an iterator to the given snarl in row_header
    std::unordered_map<stoat_vcf::Edge_t, size_t>::const_iterator find_edge(const stoat_vcf::Edge_t& edge) const;

    // Return row_header_end, an iterator to the end of row_header
    std::unordered_map<stoat_vcf::Edge_t, size_t>::const_iterator get_end_dict() const;

    // Reset row_header_end to be the end of row_header
    void set_end_dict();

    // Retrieve the index of `edge` if it exists in edge_index_dict. Otherwise, add it and return the new index.
    size_t getOrAddIndex(std::unordered_map<stoat_vcf::Edge_t, size_t>& edge_index_dict, const Edge_t& key, const size_t& size_edge_index_dict);


private:
    size_t cols_;
    size_t MaxElement;
    std::vector<uint8_t> matrix_1D;
    std::unordered_map<stoat_vcf::Edge_t, size_t> row_header;
    std::unordered_map<stoat_vcf::Edge_t, size_t>::iterator row_header_end;
};

} // end namespace stoat_vcf

#endif
