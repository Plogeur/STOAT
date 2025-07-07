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

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;

namespace stoat_vcf {

std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);
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

// Given a path, return its sample name
std::string get_sample_name_from_path(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path);

// This stores a sample name and haplotype number
struct sample_hap_t {
    std::string sample;
    std::size_t haplotype;

    const inline bool operator==(const sample_hap_t& other) const {
        return (sample==other.sample && haplotype==other.haplotype);
    }
    const inline bool operator<(const sample_hap_t& other) const {
        if (sample == other.sample) {
            return haplotype < other.haplotype;
        } else {
            return sample < other.sample;
        }
    }
};

inline std::ostream& operator<<(std::ostream& out, const sample_hap_t& sample) {
    return out << sample.sample << "#" << sample.haplotype;
}

// Given a path, return its sample name and haplotype 
sample_hap_t get_sample_and_haplotype(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path);

// equality within a given epsilon
template<typename T>
bool is_equal(T a, T b, T e = std::numeric_limits<T>::epsilon()) {
    return std::fabs(a-b) <= e;
};

enum phenotype_type_t { BINARY = 1, QUANTITATIVE, EQTL };

} // namespace stoat_vcf

#endif
