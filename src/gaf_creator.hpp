#ifndef GAF_CREATOR_HPP
#define GAF_CREATOR_HPP

#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <regex>
#include <fstream>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>

using namespace std;
using namespace bdsg;

std::pair<double, double> calcul_proportion_signi(int number_ind_group0, int number_ind_group1, double p_value);
std::string addSuffixToFilename(const std::string& filename, const std::string& suffix);
string add_suffix_to_filename(const string& filename, const string& suffix);
void write_gaf_lines(const string& sequence_name, const string& path, int length, double prop, ofstream& outfile);

vector<int> decompose_snarl(const string& snarl);
int calcul_path_length(PackedGraph& pg, const string& snarl);
void write_gaf_lines(const string& sequence_name, const string& path, int length, double prop, ofstream& outfile);
void gaf_creation(const string& input_file, std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>>& snarl_chr,
    PackedGraph& pg, const string& output_file);

#endif // GAF_CREATOR_HPP