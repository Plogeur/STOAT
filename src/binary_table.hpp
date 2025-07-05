#ifndef binary_table_HPP
#define binary_table_HPP

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "matrix.hpp"
#include "snarl_analyser.hpp"
#include "utils.hpp"
#include "stats_test.hpp"

using namespace std;

namespace stoat_vcf {

// ------------------------ Binary table ------------------------

// Write a std::string of: g0[0]:g1[1],g0[1]:g1[1],g0[2]:g1[2]...
std::string format_group_paths(const std::vector<size_t>& g0, const std::vector<size_t>& g1);

// Given two vectors of genotypes representing two groups, fill in the p-values, etc by running the relevant tests
void binary_stat_test(const std::vector<size_t>& g0, const std::vector<size_t>& g1, 
    std::string& fastfisher_p_value, std::string& chi2_p_value, std::string& group_paths,
    std::string& allele_number_str, std::string& min_row_index_str, std::string& numb_colum_str, 
    std::string& inter_group_str, std::string& average_str);

// Given two vectors of genotypes representing two groups (with length number_paths), fill them in with counts of the number of times each path is seen  
// g0 and g1 can be used in binary_stat_test()
size_t create_binary_table(
    std::vector<size_t>& g0, std::vector<size_t>& g1,
    const std::vector<bool>& binary_phenotype, 
    const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl, 
    const size_t& number_paths,
    const size_t& number_samples,
    const stoat_vcf::EdgeBySampleMatrix& matrix);

} // namespace stoat_vcf

#endif
