#ifndef QUANTITATIVE_ANALYSIS_HPP
#define QUANTITATIVE_ANALYSIS_HPP

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
#include "snarl_parser.hpp"
#include "utils.hpp"

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

// Given the number of samples (length_sample), the paths through the snarl (column_headers), the binary or quantitative phenotype of each sample (phenotype)
// and a matrix of edges in each sample
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - phenotype_filtered: the phenotypes for each genotype 
// - allele_number: the total number of alleles seen (sum of allele_paths)
// - allele_paths: the number of samples that take each path through the snarl (per column) 
// TODO: This is very similar to create_eqtl_table()
template <typename T>
std::tuple<std::vector<std::vector<double>>, std::vector<T>, size_t, std::vector<size_t>> create_quantitative_table(
    const size_t& sampleCount,
    const std::vector<std::string>& columnHeaders,
    const std::vector<T>& phenotype,
    Matrix& matrix);

// Given the number of samples (length_sample), the paths through the snarl (column_headers), and a matrix of edges in each sample,
// Return a tuple of 
// - genotypes_filtered: a matrix where each row is a sample, each column is an allele (from column_headers), counts divided by the sum of each row
// - index_used: row (samples) indices that were filled in
// - allele_number: the total number of alleles seen (sum of allele_paths)
// - allele_paths: the number of samples that take each path through the snarl (per column) 
// TODO: This is very similar to create_quantitative()
std::tuple<std::vector<std::vector<double>>, std::unordered_set<size_t>, size_t, std::vector<size_t>> create_eqtl_table(
    const size_t& length_sample,
    const std::vector<std::string>& column_headers,
    Matrix& matrix);

#endif
