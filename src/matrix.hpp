#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

class Matrix {
public:
    Matrix(uint64_t rows, uint64_t cols);
    ~Matrix()=default;
    bool operator()(uint64_t row, uint64_t col) const;
    void set(uint64_t row, uint64_t col);
    const std::vector<uint8_t>& get_matrix() const;
    const std::unordered_map<std::string, uint64_t>& get_row_header() const;
    uint64_t getMaxElement() const;
    void expandMatrix();
    void set_row_header(const std::unordered_map<std::string, uint64_t>& row_header);
    void shrink(uint64_t current_rows);

private:
    uint64_t cols_;
    uint64_t MaxElement;
    std::vector<uint8_t> matrix_1D;
    std::unordered_map<std::string, uint64_t> row_header;
};

#endif
