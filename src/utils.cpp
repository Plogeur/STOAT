#include "utils.hpp"

std::string set_precision(long double value) {
    std::ostringstream oss;
    if (value < 0.0001) {
        oss << std::scientific << std::setprecision(4) << value; // Scientific notation with 4 decimals
    } else {
        oss << std::fixed << std::setprecision(4) << value; // Fixed-point notation with 4 decimals
    }
    return oss.str();
}
