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
#include <chrono>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "matrix.hpp"

using namespace std;

// SnarlParser class declaration
class SnarlParser {
public:
    const std::string& filename;
    std::ifstream file;
    Matrix matrix;
    std::vector<std::string> sampleNames;
    
    SnarlParser(const std::string& vcf_path);
    void pushMatrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn);
    void parse_vcf(const char* vcf_filename);
    std::vector<std::string> parseHeader();
    void binary_table(const unordered_map<string, tuple<vector<string>, string, string, vector<string>>>& snarls,
                        const std::unordered_map<std::string, bool>& binary_groups,
                        const std::string& output = "output/binary_output.tsv");
    void quantitative_table(const unordered_map<string, tuple<vector<string>, string, string, vector<string>>>& snarls,
                                const std::unordered_map<std::string, double>& quantitative,
                                const std::string& output = "output/quantitative_output.tsv");
};

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
