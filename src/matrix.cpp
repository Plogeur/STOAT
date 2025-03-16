#include "matrix.hpp"

// Constructor implementation
Matrix::Matrix(size_t rows, size_t cols) 
    : cols_(cols), MaxElement(0)
{
    size_t length_matrix = (rows * cols + 7) / 8;
    row_header.rehash(rows);
    matrix_1D.reserve(length_matrix);  // Reserve capacity to avoid frequent reallocations
    matrix_1D.resize(length_matrix, 0); // Initialize with zeros
}

// Getter for matrix
const std::vector<uint8_t>& Matrix::get_matrix() const {
    return matrix_1D;
}

// Getter for row header
const std::unordered_map<std::string, size_t>& Matrix::get_row_header() const {
    return row_header;
}

// Getter row number
size_t Matrix::getMaxElement() const {
    return MaxElement;  // Convert bits back to rows
}

// Setter for row header
void Matrix::set_row_header(const std::unordered_map<std::string, size_t>& new_row_header) {
    row_header = std::move(new_row_header);
}

void Matrix::expandMatrix() {
    MaxElement *= 2;  // Double the number of elements in the matrix
    size_t new_length = matrix_1D.size() * 2;
    matrix_1D.reserve(new_length);
    matrix_1D.resize(new_length, 0); // Initialize new memory with zeros
}

// Overloaded operator() to access elements as matrix(row, col)
bool Matrix::operator()(size_t row, size_t col) const {
    size_t bitIndex = row * cols_ + col;
    size_t byteIndex = bitIndex / 8;
    size_t bitPosition = bitIndex % 8;
    if (byteIndex >= matrix_1D.size()) return false; // Bounds check to avoid out-of-range access
    return (matrix_1D[byteIndex] >> bitPosition) & 1U;
}

// Function to set a specific element (row, col) to true
void Matrix::set(size_t row, size_t col) {
    size_t bitIndex = row * cols_ + col;
    size_t byteIndex = bitIndex / 8;
    size_t bitPosition = bitIndex % 8;
    if (byteIndex >= matrix_1D.size()) return; // Bounds check to avoid out-of-range access
    matrix_1D[byteIndex] |= (1U << bitPosition);
    MaxElement = std::max(MaxElement, row + 1);  // Update MaxElement
}

void Matrix::shrink(size_t current_rows) {
    size_t new_bits = current_rows * cols_;
    size_t new_bytes = (new_bits + 7) / 8; // Compute required bytes (round up)
    matrix_1D.resize(new_bytes); // Resize
    matrix_1D.shrink_to_fit(); // Free unused capacity
}
