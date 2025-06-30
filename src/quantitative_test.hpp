#ifndef quantitative_test_HPP
#define quantitative_test_HPP

#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <iomanip>
#include <Eigen/Dense>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>  // For t-distribution
#include <boost/math/distributions/chi_squared.hpp>

#include "matrix.hpp"
#include "snarl_analyser.hpp"
#include "utils.hpp"

using namespace std;

namespace stoat_vcf{

// Linear regression function that returns a tuple of p_value, standard error (se), and beta
void linear_regression(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str);

void glm_quantitative(
    const std::vector<std::vector<double>>& df,
    const std::vector<double>& quantitative_phenotype,
    const std::vector<std::vector<double>>& covar,
    std::string& p_value_str, std::string& beta_str, 
    std::string& se_str, std::string& r2_str);

std::tuple<std::vector<std::vector<double>>, size_t, std::unordered_set<size_t>, bool, std::vector<size_t>>
process_table_quantitative(
    const size_t& number_samples,
    const std::vector<Path_traversal_t>& column_headers,
    const EdgeBySampleMatrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), the binary or quantitative phenotype of each sample (phenotype)
// and a matrix of edges in each sample
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - phenotype_filtered: the phenotypes for each genotype 
// - allele_number: the total number of alleles seen (sum of allele_paths)
// - allele_paths: the number of samples that take each path through the snarl (per column) 
template <typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>, size_t, std::vector<size_t>>
create_quantitative_table(
    const size_t& number_samples,
    const std::vector<Path_traversal_t>& column_headers,
    const std::vector<T>& phenotype,
    const EdgeBySampleMatrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), and a matrix of edges in each sample,
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - index_used: row (samples) indices that were filled in
// - allele_number: the total number of alleles seen (sum of allele_paths)
// - allele_paths: the number of samples that take each path through the snarl (per column) 
std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, size_t, std::vector<size_t>>
create_eqtl_table(
    const size_t& number_samples,
    const std::vector<Path_traversal_t>& column_headers,
    const EdgeBySampleMatrix& matrix);

}
#endif
