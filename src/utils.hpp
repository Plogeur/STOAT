#ifndef UTILS_HPP
#define UTILS_HPP

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <unordered_set>
#include <iomanip>
#include <Eigen/Dense>
#include <fstream>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "snarl_data_t.hpp"

using namespace std;
namespace stoat_vcf {

    std::string set_precision(const double& value);
    std::string set_precision_float_50(const boost::multiprecision::cpp_dec_float_50& value);

    bool is_na(const std::string& s);
    double string_to_pvalue(const std::string& p1);

    void writeSignificantTableToTSV(
        const std::vector<std::vector<double>>& table,
        const std::vector<std::string>& list_snarl,
        const std::vector<std::string>& list_samples,
        const std::string& filename);

    bool isPValueSignificant(const double& pvalue_threshold, const std::string& pvalue_str);
    void retain_indices(std::vector<double>& vec, const std::unordered_set<size_t>& indices_to_keep);
    std::vector<double> adjusted_holm(const std::vector<double>& p_values);

    template <typename T>
    std::string vectorToString(const std::vector<T>& vec);

    template <typename T>
    std::vector<T> stringToVector(const std::string& str);

} // end namespace stoat_vcf

#endif
