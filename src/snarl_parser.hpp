#pragma once

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

using namespace std;

// SnarlParser class declaration
class SnarlParser {
public:
    std::vector<std::string> sampleNames;
    Matrix matrix;

    SnarlParser(const vector<string>& sample_names, size_t num_paths_chr);
    void push_matrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn);
    void binary_table(const std::vector<std::tuple<string, vector<string>, string, string, vector<string>>>& snarls,
                        const std::unordered_map<std::string, bool>& binary_groups,
                        const std::string& output = "output/binary_gwas.tsv");
    void quantitative_table(const std::vector<std::tuple<string, vector<string>, string, string, vector<string>>>& snarls,
                                const std::unordered_map<std::string, double>& quantitative,
                                const std::string& output = "output/quantitative_gwas.tsv");
};

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

void process_vcf_batch(std::vector<bcf1_t*> records, bcf_hdr_t* hdr, SnarlParser& snarl_parser,
    std::unordered_map<std::string, size_t>& row_header_dict, std::mutex& mutex);
SnarlParser make_matrix(const std::string& vcf_filename, htsFile *vcf_file, bcf_hdr_t *hdr, bcf1_t *rec, const vector<string>& sample_names, string &chr, size_t &num_paths_ch);

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
unsigned long long int getOrAddIndex(std::unordered_map<std::string, unsigned long long int>& orderedMap, const std::string& key, unsigned long long int lengthOrderedMap);

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s);

// Function to split a string by a delimiter
std::vector<std::string> split(const std::string& str, char delimiter);

std::vector<int> extractGenotype(const std::string& genotypeStr);

std::vector<std::string> extractATField(const std::string& infoField);

std::vector<std::string> split(const std::string& s, char delimiter);

// Function to determine and extract an integer from the string
std::pair<int, std::string> determine_str(const std::string& s, size_t length_s, size_t i);

// Function to decompose a list of snarl strings
const std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst);

std::vector<int> identify_correct_path(const std::vector<std::string>& decomposed_snarl, 
                                        const std::unordered_map<std::string, 
                                        size_t>& row_headers_dict, 
                                        const Matrix& matrix, 
                                        const size_t num_cols);
