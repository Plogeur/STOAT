#include "matrix.hpp"

// Constructor implementation
Matrix::Matrix(uint64_t rows, uint64_t cols) 
    : cols_(cols) {

    uint64_t length_matrix = (rows * cols + 7) / 8;
    MaxElement = (length_matrix * 8) / cols_; // get the number of element in the matrix
    row_header.rehash(rows);
    matrix_1D.reserve(length_matrix);  // Reserve capacity to avoid frequent reallocations
    matrix_1D.resize(length_matrix, 0); // Initialize with zeros
}

// Getter for matrix
const std::vector<uint8_t>& Matrix::get_matrix() const {
    return matrix_1D;
}

// Getter for row header
const std::unordered_map<std::string, uint64_t>& Matrix::get_row_header() const {
    return row_header;
}

// Getter row number
uint64_t Matrix::getMaxElement() const {
    return MaxElement;  // Convert bits back to rows
}

// Setter for row header
void Matrix::set_row_header(const std::unordered_map<std::string, uint64_t>& new_row_header) {
    row_header = std::move(new_row_header);
}

void Matrix::expandMatrix() {
    MaxElement *= 2;  // Double the number of elements in the matrix
    uint64_t new_length = matrix_1D.size() * 2;
    matrix_1D.reserve(new_length);
    matrix_1D.resize(new_length, 0); // Initialize new memory with zeros
}

// Overloaded operator() to access elements as matrix(row, col)
bool Matrix::operator()(uint64_t row, uint64_t col) const {
    uint64_t bitIndex = row * cols_ + col;
    uint64_t byteIndex = bitIndex / 8;
    uint64_t bitPosition = bitIndex % 8;
    if (byteIndex >= matrix_1D.size()) return false; // Bounds check to avoid out-of-range access
    return (matrix_1D[byteIndex] >> bitPosition) & 1U;
}

// Function to set a specific element (row, col) to true
void Matrix::set(uint64_t row, uint64_t col) {
    uint64_t bitIndex = row * cols_ + col;
    uint64_t byteIndex = bitIndex / 8;
    uint64_t bitPosition = bitIndex % 8;
    if (byteIndex >= matrix_1D.size()) return; // Bounds check to avoid out-of-range access
    matrix_1D[byteIndex] |= (1U << bitPosition);
}

void Matrix::shrink(uint64_t current_rows) {
    uint64_t new_bits = current_rows * cols_;
    uint64_t new_bytes = (new_bits + 7) / 8; // Compute required bytes (round up)
    matrix_1D.resize(new_bytes); // Resize
    matrix_1D.shrink_to_fit(); // Free unused capacity
}
