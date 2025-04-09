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

// Struct to store parsed eQTL data// QTL struct
struct QTL {
    std::vector<std::string> sample_ids;
    std::vector<std::string> gene_ids;
    std::vector<std::vector<double>> expression_matrix;

    // Default constructor
    QTL() = default;
    
    QTL(const std::vector<std::string>& sample_ids,
        const std::vector<std::string>& gene_ids,
        const std::vector<std::vector<double>>& expression_matrix)
        : sample_ids(sample_ids), gene_ids(gene_ids), expression_matrix(expression_matrix) {}
};

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
};

KinshipMatrix parseKinshipMatrix(const std::string& filename);
QTL parseExpressionFile(const std::string& filename);

void check_format_covariate(const std::string& filename);

template <typename T>
void check_phenotype_covariate(const std::unordered_map<std::string, T>& phenotype, 
    const std::unordered_map<std::string, std::vector<double>>& covariates);

std::unordered_map<std::string, std::vector<double>> parseCovariate(const std::string& filename);

// Parses the group file and fills the group_0 and group_1 maps with sample data.
std::unordered_map<std::string, bool> parse_binary_pheno(const std::string& binary_pheno);

// Parses the phenotype file and returns a map with IID as keys and PHENO as float values.
std::unordered_map<std::string, double> parse_quantitative_pheno(const std::string& qunatitative_pheno);

std::tuple<std::vector<std::string>, htsFile*, bcf_hdr_t*, bcf1_t*> parseHeader(const std::string& file_path);

std::unordered_set<std::string> parse_chromosome_reference(const string& file_path);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys);

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> parse_snarl_path(const std::string& path_file);

void check_format_quantitative_phenotype(const std::string& file_path);
void check_format_binary_phenotype(const std::string& file_path);
void check_file(const std::string& file_path);

#endif
