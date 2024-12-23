#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

class Matrix {
public:
    Matrix(size_t rows, size_t cols);
    bool operator()(size_t row, size_t col) const;
    void set(size_t row, size_t col);
    const std::vector<uint8_t>& get_matrix() const;
    const std::unordered_map<std::string, size_t>& get_row_header() const;
    size_t getRows() const;
    void expandMatrix();
    void set_row_header(const std::unordered_map<std::string, size_t>& row_header);

private:
    size_t rows_;
    size_t cols_;
    size_t default_row_number;
    std::vector<uint8_t> matrix_1D;
    std::unordered_map<std::string, size_t> row_header;
};

#endif
