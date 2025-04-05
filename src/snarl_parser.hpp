#ifndef SNARL_PARSER_HPP
#define SNARL_PARSER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <thread>
#include <mutex>
#include <chrono>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "matrix.hpp"
#include "arg_parser.hpp"

using namespace std;

// SnarlParser class declaration
class SnarlParser {
public:
    std::vector<std::string> sampleNames;
    Matrix matrix;

    SnarlParser(const vector<string>& sample_names, size_t num_paths_chr);
    ~SnarlParser()=default;
    void push_matrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn);
    
    void binary_table(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls,
                        const std::unordered_map<std::string, bool>& binary_groups, const string &chr,
                        const std::unordered_map<std::string, std::vector<double>>& covar,
                        const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
                        const KinshipMatrix& kinship, std::ofstream& outf);

    void quantitative_table(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls,
                            const std::unordered_map<std::string, double>& quantitative_phenotype, const string &chr,
                            const std::unordered_map<std::string, std::vector<double>>& covar,
                            const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector,
                            const KinshipMatrix& kinship, std::ofstream& outf);

    void create_bim_bed(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls, 
                                    string chromosome, const std::string& output_bim, const std::string& output_bed);
    std::vector<int> create_table_short_path(const std::string& list_path_snarl);

};

// void chromosome_chuck_eqtl(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
//     const std::vector<std::string> &list_samples,
//     unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
//     const vector<QTLRecord> pheno, std::ofstream& outf);

bool check_MAF_threshold_binary(const std::vector<std::vector<int>>& df, const double& maf);
bool check_MAF_threshold_quantitative(const std::unordered_map<std::string, std::vector<int>>& df, const double& maf);

void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples, 
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    const unordered_map<string, bool>& pheno, std::unordered_map<std::string, std::vector<double>> covar, 
    const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
    const KinshipMatrix& kinship, const std::string& output_binary);

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    const unordered_map<string, double>& pheno, std::unordered_map<std::string, std::vector<double>> covar,
    const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
    const KinshipMatrix& kinship, const std::string& output_quantitive);

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    const unordered_map<string, double>& pheno, string output_dir);

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

void create_fam(const std::vector<std::pair<std::string, int>> &pheno, 
    const std::string& output_path);

std::tuple<SnarlParser, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const vector<string>& sample_names, string &chr, size_t &num_paths_ch);

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<std::string, size_t>& orderedMap, const std::string& key, size_t lengthOrderedMap);

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s);

// Function to determine and extract an integer from the string
std::pair<int, std::string> determine_str(const std::string& s, size_t length_s, size_t i);

// Function to decompose a list of snarl strings
const std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst);

std::vector<int> identify_correct_path(const std::vector<std::string>& decomposed_snarl, 
                                        const std::unordered_map<std::string, 
                                        size_t>& row_headers_dict, 
                                        const Matrix& matrix, 
                                        const size_t num_cols);

            std::unordered_map<std::string, std::vector<double>> convertBinaryGroups(
                const std::unordered_map<std::string, bool>& binary_groups);

#endif