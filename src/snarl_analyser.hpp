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

#include "arg_parser.hpp"
#include "matrix.hpp"
#include "snarl_data_t.hpp"
#include "binary_table.hpp"
#include "quantitative_table.hpp"
#include "utils.hpp"

using namespace std;

namespace stoat_vcf {

/// A SnarlAnalyser stores a bit matrix of samples (columns) vs edges they take (rows), taken from a VCF
/// It also has a vector of sample names 
class SnarlAnalyser {
public:
    // caca change to private
    std::vector<std::string> sampleNames;
    EdgeBySampleMatrix matrix;

    SnarlAnalyser(const std::vector<std::string>& sample_names, size_t num_paths_chr);
    ~SnarlAnalyser()=default;
    void push_matrix(const Edge_t& EdgePath, std::unordered_map<stoat_vcf::Edge_t, size_t>& edge_dict, size_t indexColumn);
    
    void binary_table(const std::vector<Snarl_data_t>& snarls,
        const std::vector<bool>& binary_phenotype, const std::string& chr,
        const std::vector<std::vector<double>>& covar,
        const double& maf, const size_t& num_threads, 
        const double& table_threshold, const std::string& output_dir, std::ofstream& outf);

    /// Similar to binary_table, get the genotypes and write the tsv output
    void quantitative_table(const std::vector<Snarl_data_t>& snarls,
                            const std::vector<double>& quantitative_phenotype, const std::string &chr,
                            const std::vector<std::vector<double>>& covar,
                            const double& maf, const size_t& num_threads, 
                            const double& table_threshold, const std::string& output_dir, std::ofstream& outf);

    /// For each snarl, write a bim file and a bed file (PLINK formats)
    void create_bim_bed(const std::vector<Snarl_data_t>& snarls, 
        std::string chromosome, std::ofstream& outbim, std::ofstream& outbed);

    /// Similar to binary_table and quantitative_table, get the genotype and write the tsv output
    void eqtl_table(const std::vector<Snarl_data_t>& snarls,
        const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& eqtl,
        const std::string& chr, const std::vector<std::vector<double>>& covar,
        const double& maf, const size_t& num_threads, 
        const double& table_threshold, const std::string& regression_dir,
        const size_t& windows_gene_threshold, std::ofstream& outf);

    /// Given a list of paths in a snarl, return vectors of counts of each sample taking an allele for the two alleles with the highest counts over all samples
    /// Each vector returned corresponds to one allele, each entry in the vector is a sample, the value is a count of the number of times a sample takes the path/allele (probably binary?)
    std::pair<std::vector<size_t>, std::vector<size_t>> create_table_short_path(const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl);
};

/// Return true if any column exceeds the MAF threshold 
bool check_MAF_threshold_quantitative(const std::vector<std::vector<double>>& df, const double& maf);

/// Go through the vcf by chromosome, parse it to get a matrix of genotypes (SnarlParser of edges), then write the binary table
void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples, 
    const std::unordered_map<std::string, std::vector<Snarl_data_t>> &snarl_chr,
    const std::vector<bool>& pheno, std::vector<std::vector<double>> covar, 
    const double& maf,
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const std::string& output_binary);

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>> &snarl_chr,
    const std::vector<double>& pheno, std::vector<std::vector<double>> covar,
    const double& maf,
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const std::string& output_quantitative);

void chromosome_chuck_eqtl(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& snarl_chr,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>>& eqtl_map,
    const std::vector<std::vector<double>>& covar,
    const double& maf,
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const size_t& windows_gene_threshold, 
    const std::string& out_eqtl);

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& snarl_chr,
    const std::string& output_dir);

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

std::vector<size_t> found_gene_snarl(
    const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold);

void create_fam(const std::vector<std::pair<std::string, int>> &pheno, 
    const std::string& output_path);

/// Make a SnarlParser representing the genotypes in a vcf and the pointers to the vcf but advanced to the end of the chromosome?
std::tuple<SnarlAnalyser, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const std::vector<std::string>& sample_names, std::string &chr, size_t &num_paths_ch);

// Retrieve the index of `edge` if it exists in edge_index_dict. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<stoat_vcf::Edge_t, size_t>& edge_index_dict, const Edge_t& key, const size_t& size_edge_index_dict);

// Function to determine and extract an node id from the std::string
inline size_t extract_node_id(const std::string& s, size_t length_s, size_t& i);

// Decompose path Path_traversal_t to vector Edge_t
std::vector<stoat_vcf::Edge_t> decompose_path_to_edges(const Path_traversal_t& s);

// Decompose a list of paths Path_traversal_t into a vector of Edge_t
const std::vector<std::vector<stoat_vcf::Edge_t>> decompose_path_list_path(const std::vector<stoat_vcf::Path_traversal_t>& list_paths);

// Decompose a list of paths std::string into a vector of Edge_t
const std::vector<std::vector<stoat_vcf::Edge_t>> decompose_path_list_str(const std::vector<std::string>& list_paths);

// Decompose path std::string to vector Edge_t
std::vector<stoat_vcf::Edge_t> decompose_path_str_to_edge(const std::string& s);

/// Given a path through the snarl, a matrix of edges for each sample/haplotype, and the number of columns (samples/haplotypes),
/// return the columns for which all edges (rows) in the path are set, i.e. the haplotypes with the given path.
std::vector<size_t> identify_path(
    const std::vector<stoat_vcf::Edge_t>& list_edge_path,
    const EdgeBySampleMatrix& matrix,
    const size_t num_cols);

std::vector<std::vector<size_t>> transpose_matrix(const std::vector<std::vector<size_t>>& matrix);

/// Set major_index_1 and major_index2 to be the indices of the largest and second largest values in vec
void find_two_largest_indices(const std::vector<size_t>& vec, size_t& major_index_1, size_t& major_index_2);

} //end stoat_vcf namespace

#endif
