#include "matrix.hpp"

namespace stoat_vcf {

// Constructor implementation
EdgeBySampleMatrix::EdgeBySampleMatrix(const std::vector<std::string>& sampleNames, size_t rows, size_t cols) : cols_(cols), sampleNames(sampleNames) {

    size_t length_matrix = (rows * cols + 7) / 8;
    MaxElement = (length_matrix * 8) / cols_; // get the number of element in the matrix
    row_header.rehash(rows);
    matrix_1D.reserve(length_matrix); // Reserve capacity to avoid frequent reallocations
    matrix_1D.resize(length_matrix, 0); // Initialize with zeros
}

// Getter for row header
std::unordered_map<stoat_vcf::Edge_t, size_t>::const_iterator EdgeBySampleMatrix::find_edge(const stoat_vcf::Edge_t& edge_s) const {
    return row_header.find(edge_s);
}

// Getter for row header
std::unordered_map<stoat_vcf::Edge_t, size_t>::const_iterator EdgeBySampleMatrix::get_end_dict() const {
    return row_header_end;
}

// Retrieve the index of `key` if it exists in the dict. Otherwise, add it and return the new index.
size_t EdgeBySampleMatrix::getOrAddIndex(std::unordered_map<stoat_vcf::Edge_t, size_t>& edge_index_dict, const Edge_t& key, const size_t& size_edge_index_dict) {
    auto it = edge_index_dict.find(key);
    if (it != edge_index_dict.end()) {
        return it->second;
    } else {
        size_t newIndex = size_edge_index_dict;
        edge_index_dict[key] = newIndex;
        return newIndex;
    }
}


// Add True to the matrix if edge is found
void EdgeBySampleMatrix::push_matrix(const Edge_t& EdgePath, std::unordered_map<stoat_vcf::Edge_t, size_t>& edge_index_dict, size_t indexColumn) {

    size_t lengthOrderedMap = edge_index_dict.size();
    size_t idxSnarl = getOrAddIndex(edge_index_dict, EdgePath, lengthOrderedMap);
    size_t currentRowsNumber = getMaxElement();

    if (lengthOrderedMap > currentRowsNumber - 1) {
        expandMatrix();
    }

    set(idxSnarl, indexColumn);
}


// Getter for row header
void EdgeBySampleMatrix::set_end_dict() {
    row_header_end = row_header.end();
}

// Getter row number
size_t EdgeBySampleMatrix::getMaxElement() const {
    return MaxElement;  // Convert bits back to rows
}

// Setter for row header
// TODO: This could be done when things are added to the matrix in push_matrix
void EdgeBySampleMatrix::set_row_header(const std::unordered_map<stoat_vcf::Edge_t, size_t>& new_row_header) {
    row_header = std::move(new_row_header);
}

void EdgeBySampleMatrix::expandMatrix() {
    MaxElement *= 2;  // Double the number of elements in the matrix
    size_t new_length = matrix_1D.size() * 2;
    matrix_1D.reserve(new_length);
    matrix_1D.resize(new_length, 0); // Initialize new memory with zeros
}

// Overloaded operator() to access elements as matrix(row, col)
bool EdgeBySampleMatrix::operator()(size_t row, size_t col) const {
    size_t bitIndex = row * cols_ + col;
    size_t byteIndex = bitIndex / 8;
    size_t bitPosition = bitIndex % 8;
    // Bounds check to avoid out-of-range access
    // if (byteIndex >= matrix_1D.size()) return false;
    return (matrix_1D[byteIndex] >> bitPosition) & 1U;
}

// Function to set a specific element (row, col) to true
void EdgeBySampleMatrix::set(size_t row, size_t col) {
    size_t bitIndex = row * cols_ + col;
    size_t byteIndex = bitIndex / 8;
    size_t bitPosition = bitIndex % 8;
    // Bounds check to avoid out-of-range access
    // if (byteIndex >= matrix_1D.size()) return;
    matrix_1D[byteIndex] |= (1U << bitPosition);
}

void EdgeBySampleMatrix::shrink(size_t current_rows) {
    size_t new_bits = current_rows * cols_;
    size_t new_bytes = (new_bits + 7) / 8; // Compute required bytes (round up)
    matrix_1D.resize(new_bytes); // Resize
    matrix_1D.shrink_to_fit(); // Free unused capacity
}

} // end namespace stoat_vcf
