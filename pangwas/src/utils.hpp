#ifndef PANGWAS_UTILS_HPP_INCLUDED
#define PANGWAS_UTILS_HPP_INCLUDED

#include <string>
#include <cmath>
#include <handlegraph/path_position_handle_graph.hpp>

using namespace std;

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

#endif
