#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>


// A class to store a 2d bit-matrix
// Rows represent edges and the index of each edge can be found from the row_header
// Columns represent samples/haplotypes
class Matrix {
public:
    Matrix(size_t rows, size_t cols);
    ~Matrix()=default;

    // Operator to get the value
    bool operator()(size_t row, size_t col) const;

    // Set this value to true
    void set(size_t row, size_t col);

    // Get the matrix itself
    const std::vector<uint8_t>& get_matrix() const;

    // Get the row_header
    const std::unordered_map<std::string, size_t>& get_row_header() const;

    // Get the maximum index into the vector representing the matrix
    size_t getMaxElement() const;
    
    // Double the size of the matrix
    void expandMatrix();

    // Reset row_header
    void set_row_header(const std::unordered_map<std::string, size_t>& row_header);

    // Shrink to use the minimum amount of memory possible allowing current_rows
    void shrink(size_t current_rows);

    // Return an iterator to the given snarl in row_header
    // TODO: I think this should be edge not snarl
    std::unordered_map<std::string, size_t>::const_iterator find_snarl(const std::string& snarl) const;

    // Return row_header_end, an iterator to the end of row_header
    std::unordered_map<std::string, size_t>::const_iterator get_end_dict() const;

    // Reset row_header_end to be the end of row_header
    void set_end_dict();

private:
    size_t cols_;
    size_t MaxElement;
    std::vector<uint8_t> matrix_1D;
    std::unordered_map<std::string, size_t> row_header;
    std::unordered_map<std::string, size_t>::iterator row_header_end;
};

#endif
