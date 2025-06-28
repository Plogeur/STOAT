#include "matrix.hpp"

// Constructor implementation
EdgeBySampleMatrix::EdgeBySampleMatrix(size_t rows, size_t cols) : cols_(cols) {

    size_t length_matrix = (rows * cols + 7) / 8;
    MaxElement = (length_matrix * 8) / cols_; // get the number of element in the matrix
    row_header.rehash(rows);
    matrix_1D.reserve(length_matrix); // Reserve capacity to avoid frequent reallocations
    matrix_1D.resize(length_matrix, 0); // Initialize with zeros
}

// Getter for matrix
const std::vector<uint8_t>& EdgeBySampleMatrix::get_matrix() const {
    return matrix_1D;
}

// Getter for row header
std::unordered_map<Edge_t, size_t>::const_iterator EdgeBySampleMatrix::find_edge(const &Edge_t edge) const {
    return row_header.find(edge);
}

// Getter for row header
std::unordered_map<Edge_t, size_t>::const_iterator EdgeBySampleMatrix::get_end_dict() const {
    return row_header_end;
}

// Getter for row header
void EdgeBySampleMatrix::set_end_dict() {
    row_header_end = row_header.end();
}

// Getter for row header
const std::unordered_map<Edge_t, size_t>& EdgeBySampleMatrix::get_row_header() const {
    return row_header;
}

// Getter row number
size_t EdgeBySampleMatrix::getMaxElement() const {
    return MaxElement;  // Convert bits back to rows
}

// Setter for row header
void EdgeBySampleMatrix::set_row_header(const std::unordered_map<Edge_t, size_t>& new_row_header) {
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
