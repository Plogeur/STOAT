#ifndef ARG_PARSER_HPP
#define ARG_PARSER_HPP

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <map>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <stdexcept>
#include <regex>
#include <Eigen/Dense>

#include <htslib/vcf.h>
#include <htslib/hts.h>

using namespace std;

// KinshipMatrix struct
struct KinshipMatrix {
    std::vector<std::string> ids;
    std::vector<std::vector<double>> matrix;

    // Default constructor
    KinshipMatrix() = default;

    // Parameterized constructor
    KinshipMatrix(const std::vector<std::string>& ids,
                  const std::vector<std::vector<double>>& matrix)
        : ids(ids), matrix(matrix) {}
    
    const bool empty() const;
};

KinshipMatrix parseKinshipMatrix(const std::string& filename);
std::unordered_map<std::string, std::vector<double>> parse_qtl_file(const std::string& filename);
std::unordered_map<std::string, std::tuple<std::string, int, int>> parse_gene_positions(const std::string& filename);

void check_qtl_gene_position(const std::unordered_map<std::string, std::vector<double>>& qtl, const std::unordered_map<std::string, std::tuple<std::string, int, int>>& gene_position);
void check_format_covariate(const std::string& filename);

template <typename T>
void check_phenotype_covariate(const std::unordered_map<std::string, T>& phenotype, 
    const std::unordered_map<std::string, std::vector<double>>& covariates);

std::unordered_map<std::string, std::vector<double>> parse_covariates(
    const std::string& filename, const std::vector<std::string>& covar_names);

// Parses the group file and fills the group_0 and group_1 maps with sample data.
std::vector<bool> parse_binary_pheno(
    const std::string& file_path,
    const std::vector<std::string>& list_samples);

// Parses the phenotype file and returns a map with IID as keys and PHENO as float values.
std::vector<double> parse_quantitative_pheno(
    const std::string& file_path, 
    const std::vector<std::string>& list_samples);

std::tuple<std::vector<std::string>, htsFile*, bcf_hdr_t*, bcf1_t*> parseHeader(const std::string& file_path);

std::unordered_set<std::string> parse_chromosome_reference(const string& file_path);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys);

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>> parse_snarl_path(const std::string& path_file);

void check_file(const std::string& file_path);

#endif
