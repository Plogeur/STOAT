#ifndef snarl_analyser_HPP
#define snarl_analyser_HPP

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
#include <future>
#include <chrono>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "matrix.hpp"
#include "arg_parser.hpp"
#include "snarl_data_t.hpp"

using namespace std;

// SnarlAnalyser class declaration
class SnarlAnalyser {
public:
    // caca change to private
    std::vector<std::string> sampleNames;
    EdgeBySampleMatrix matrix;

    SnarlAnalyser(const vector<string>& sample_names, size_t num_paths_chr);
    ~SnarlAnalyser()=default;
    void push_matrix(const Edge_t& EdgePath, std::unordered_map<Edge_t, size_t>& edge_dict, size_t indexColumn);
    
    void binary_table(const std::vector<Snarl_data_t>& snarls,
        const std::vector<bool>& binary_phenotype, const std::string& chr,
        const std::vector<std::vector<double>>& covar,
        const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
        const double& table_threshold, const std::string& output_dir, std::ofstream& outf);

    void quantitative_table(const std::vector<Snarl_data_t>& snarls,
                            const std::vector<double>& quantitative_phenotype, const string &chr,
                            const std::vector<std::vector<double>>& covar,
                            const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
                            const double& table_threshold, const std::string& output_dir, std::ofstream& outf);

    void create_bim_bed(const std::vector<Snarl_data_t>& snarls, 
        string chromosome, std::ofstream& outbim, std::ofstream& outbed);

    void eqtl_table(
        const std::vector<Snarl_data_t>& snarls,
        const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& eqtl,
        const std::string& chr, const std::vector<std::vector<double>>& covar,
        const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
        const double& table_threshold, const std::string& regression_dir,
        const size_t& windows_gene_threshold, std::ofstream& outf);

    std::pair<std::vector<size_t>, std::vector<size_t>> create_table_short_path(const vector<Path_traversal_t>& list_path_snarl);
};

bool check_MAF_threshold_quantitative(const std::vector<std::vector<double>>& df, const double& maf);

void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples, 
    const std::unordered_map<std::string, std::vector<Snarl_data_t>> &snarl_chr,
    const std::vector<bool>& pheno, std::vector<std::vector<double>> covar, 
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const std::string& output_binary);

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>> &snarl_chr,
    const vector<double>& pheno, std::vector<std::vector<double>> covar,
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const std::string& output_quantitative);

void chromosome_chuck_eqtl(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& snarl_chr,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>>& eqtl_map,
    const std::vector<std::vector<double>>& covar,
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const size_t& windows_gene_threshold, 
    const std::string& out_eqtl);

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& snarl_chr,
    const string& output_dir);

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

std::vector<size_t> found_gene_snarl(
    const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold);

void create_fam(const std::vector<std::pair<std::string, int>> &pheno, 
    const std::string& output_path);

std::tuple<SnarlAnalyser, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const vector<string>& sample_names, string &chr, size_t &num_paths_ch);

// Retrieve the index of `key` if it exists in edge_index_dict. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<Edge_t, size_t>& edge_index_dict, const Edge_t& key, const size_t& size_edge_index_dict);

// Function to determine and extract an integer from the string
inline size_t extract_node_id(const std::string& s, size_t length_s, size_t& i);

// Decompose path Path_traversal_t to vector Edge_t
std::vector<Edge_t> decompose_path_to_edges(const Path_traversal_t& s);

// Decompose a list of paths Path_traversal_t into a vector of Edge_t
const std::vector<std::vector<Edge_t>> decompose_path_list_str(const std::vector<Path_traversal_t>& list_paths);

// Decompose path string to vector Edge_t
vector<Edge_t> decompose_path_str_to_edge(const std::string& s);

// Decompose a list of paths str into a vector of Edge_t
const std::vector<std::vector<Edge_t>> decompose_path_list_str(const std::vector<std::string>& list_paths);

std::vector<size_t> identify_path(
    const std::vector<Edge_t>& list_edge_path,
    const EdgeBySampleMatrix& matrix,
    const size_t num_cols);

std::vector<std::vector<size_t>> transpose_matrix(const std::vector<std::vector<size_t>>& matrix);

void find_two_largest_indices(const std::vector<size_t>& vec, size_t& major_index_1, size_t& major_index_2);

#endif