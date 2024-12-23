#ifndef ARG_PARSERS_HPP
#define ARG_PARSERS_HPP

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <stdexcept>

// Parses the group file and fills the group_0 and group_1 maps with sample data.
std::unordered_map<std::string, bool> parse_binary_pheno(const std::string& binary_pheno);

// Parses the phenotype file and returns a map with IID as keys and PHENO as float values.
std::unordered_map<std::string, float> parse_quantitative_pheno(const std::string& qunatitative_pheno);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys);

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<std::string>> parse_snarl_path(const std::string& path_file);

void check_format_paths_snarl(const std::string& file_path);
void check_format_quantitative_phenotype(const std::string& file_path);
void check_format_binary_phenotype(const std::string& file_path);
void check_format_vcf_file(const std::string& file_path);

std::vector<std::string> parseHeader(const std::string& file_path);

#endif
