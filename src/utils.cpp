#include "utils.hpp"

std::string set_precision(double value) {
    std::ostringstream oss;
    if (value < 0.0001) {
        oss << std::scientific << std::setprecision(4) << value; // Scientific notation with 4 decimals
    } else {
        oss << std::fixed << std::setprecision(4) << value; // Fixed-point notation with 4 decimals
    }
    return oss.str();
}

bool is_na(const std::string& s) {
    return s.empty() || s == "NA";
}

double mean_pvalue_from_strings(const std::string& p1, const std::string& p2) {
    bool na1 = is_na(p1);
    bool na2 = is_na(p2);

    if (!na1 && !na2) {
        return (std::stod(p1) + std::stod(p2)) / 2.0;
    } else if (!na1) {
        return std::stod(p1);
    } else if (!na2) {
        return std::stod(p2);
    } else {
        return 1.0;
    }
}

double string_to_pvalue(const std::string& p1) {
    bool na1 = is_na(p1);

    if (!na1) {
        return std::stod(p1);
    } else {
        return 1.0;
    }
}