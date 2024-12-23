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

// Parse file
void parse_sex(const std::string& path, std::unordered_map<std::string, int>& sex_map);
void parse_pheno(const std::string &pheno_path, std::unordered_map<std::string, int> &pheno);

void check_match_samples(const std::unordered_map<std::string, int>& map, const std::vector<std::string>& keys);

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<std::string>> parse_snarl_path(const std::string& path_file);

void check_format_paths_snarl(const std::string& file_path);
void check_format_phenotype(const std::string& file_path);
void check_format_vcf_file(const std::string& file_path);

std::vector<std::string> parseHeader(const std::string& file_path);

#endif
