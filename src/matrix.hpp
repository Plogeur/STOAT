#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdint>

#include "snarl_data_t.hpp"

class EdgeBySampleMatrix {
public:
    EdgeBySampleMatrix(size_t rows, size_t cols);
    ~EdgeBySampleMatrix()=default;
    bool operator()(size_t row, size_t col) const;
    void set(size_t row, size_t col);
    const std::vector<uint8_t>& get_matrix() const;
    const std::unordered_map<Edge_t, size_t>& get_row_header() const;
    size_t getMaxElement() const;
    void expandMatrix();
    void set_row_header(const std::unordered_map<Edge_t, size_t>& row_header);
    void shrink(size_t current_rows);
    std::unordered_map<Edge_t, size_t>::const_iterator find_edge(const Edge_t& edge) const;
    std::unordered_map<Edge_t, size_t>::const_iterator get_end_dict() const;
    void set_end_dict();

private:
    size_t cols_;
    size_t MaxElement;
    std::vector<uint8_t> matrix_1D;
    std::unordered_map<Edge_t, size_t> row_header;
    std::unordered_map<Edge_t, size_t>::iterator row_header_end;
};

#endif
