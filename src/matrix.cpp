#include "matrix.hpp"

// Constructor implementation
Matrix::Matrix(size_t rows, size_t cols) 
    : cols_(cols), default_row_number(rows)
{
    size_t length_matrix = default_row_number * cols_;
    row_header.rehash(rows); 
    matrix_1D.reserve((length_matrix + 7) / 8);  // Round up to account for any leftover bits
    matrix_1D.resize((length_matrix + 7) / 8, 0);
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
size_t Matrix::getRows() const {
    return matrix_1D.size() * 8;
}

// Setter for row header
void Matrix::set_row_header(const std::unordered_map<std::string, size_t>& row_header) {
    this->row_header = row_header;
}

void Matrix::expandMatrix() {
    size_t current_elements = matrix_1D.size() * 8;  // Number of booleans (since each uint8_t stores 8 booleans)
    size_t new_elements = (current_elements / 8) + (default_row_number * cols_); // double capacity 
    size_t new_matrix_size = (new_elements + 7) / 8;
    matrix_1D.reserve(new_matrix_size);
    matrix_1D.resize(new_matrix_size, 0);
}

// Overloaded operator() to access elements as matrix(row, col)
bool Matrix::operator()(size_t row, size_t col) const {
    size_t bitIndex = row * cols_ + col;
    size_t byteIndex = bitIndex / 8;
    size_t bitPosition = bitIndex % 8;
    return (matrix_1D[byteIndex] >> bitPosition) & 1;
}

// Function to set a specific element (row, col) modify false to true (only this ways)
void Matrix::set(size_t row, size_t col) {
    size_t bitIndex = row * cols_ + col;
    size_t byteIndex = bitIndex / 8;
    size_t bitPosition = bitIndex % 8;
    matrix_1D[byteIndex] |= (1U << bitPosition);
}

void Matrix::shrink(size_t current_rows) {

    size_t new_bits = current_rows * cols_; 
    size_t new_bytes = (new_bits + 7) / 8; // Compute required bytes (round up)

    matrix_1D.resize(new_bytes);
    matrix_1D.shrink_to_fit(); // Shrink to fit the new size
}
