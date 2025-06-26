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
#include <future>
#include <chrono>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#include "matrix.hpp"
#include "arg_parser.hpp"

using namespace std;

// SnarlParser class declaration
// A SnarlParser stores a bit matrix of samples (columns) vs edges they take (rows), taken from a VCF
// It also has a vector of sample names 
// TODO: I think this is more of a vcf parser than a snarl parser
class SnarlParser {
public:
    // caca change to private
    std::vector<std::string> sampleNames;
    Matrix matrix;

    // Create a snarl parser
    // Initialize the matrix to have num_paths_chr*4 rows and sample_names.size()*2 columns
    // The number of rows will change as more edges are seen
    SnarlParser(const vector<string>& sample_names, size_t num_paths_chr);
    ~SnarlParser()=default;

    // Given a string representing a path (edge?), the row dictionary from the Matrix, and the index of the column (sample), set the
    // matrix bit 
    void push_matrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn);
    
    // Given a vector of snarls (snarl, vector of paths, start index, end index, variant type),
    // the binary phenotype of each sample, the covariate matrix, MAF threshold, kinship matrix, number of threads,
    // a p-value threshold?, output directory, and output file,
    // Get the genotype/phenotype counts, write the tsv output 
    void binary_table(const std::vector<std::tuple<std::string, std::vector<std::string>, size_t, size_t, std::vector<std::string>>>& snarls,
        const std::vector<bool>& binary_phenotype, const std::string& chr,
        const std::vector<std::vector<double>>& covar,
        const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
        const double& table_threshold, const std::string& output_dir, std::ofstream& outf);

    void quantitative_table(const std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>& snarls,
                            const std::vector<double>& quantitative_phenotype, const string &chr,
                            const std::vector<std::vector<double>>& covar,
                            const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
                            const double& table_threshold, const std::string& output_dir, std::ofstream& outf);

    void create_bim_bed(const std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>& snarls, 
        string chromosome, std::ofstream& outbim, std::ofstream& outbed);

    void eqtl_table(
        const std::vector<std::tuple<std::string, std::vector<std::string>, size_t, size_t, std::vector<std::string>>>& snarls,
        const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& eqtl,
        const std::string& chr, const std::vector<std::vector<double>>& covar,
        const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
        const double& table_threshold, const std::string& regression_dir,
        const size_t& windows_gene_threshold, std::ofstream& outf);

    std::pair<std::vector<size_t>, std::vector<size_t>> create_table_short_path(const vector<std::string>& list_path_snarl);
};

bool check_MAF_threshold_quantitative(const std::vector<std::vector<double>>& df, const double& maf);


// Go through the vcf by chromosome, parse it to get a matrix of genotypes (SnarlParser of edges), then write the binary table
void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples, 
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>, size_t, size_t, std::vector<std::string>>>> &snarl_chr,
    const std::vector<bool>& pheno, std::vector<std::vector<double>> covar, 
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const std::string& output_binary);

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>, size_t, size_t, std::vector<std::string>>>> &snarl_chr,
    const vector<double>& pheno, std::vector<std::vector<double>> covar,
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const std::string& output_quantitative);

void chromosome_chuck_eqtl(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>, size_t, size_t, std::vector<std::string>>>>& snarl_chr,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>>& eqtl_map,
    const std::vector<std::vector<double>>& covar,
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const double& table_threshold, 
    const std::string& regression_dir, const size_t& windows_gene_threshold, 
    const std::string& out_eqtl);

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<std::string>, size_t, size_t, std::vector<std::string>>>>& snarl_chr,
    const string& output_dir);

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path);

std::vector<size_t> found_gene_snarl(
    const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold);

void create_fam(const std::vector<std::pair<std::string, int>> &pheno, 
    const std::string& output_path);

// Make a SnarlParser representing the genotypes in a vcf and the pointers to the vcf but advanced to the end of the chromosome?
std::tuple<SnarlParser, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const vector<string>& sample_names, string &chr, size_t &num_paths_ch);

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<std::string, size_t>& orderedMap, const std::string& key, size_t lengthOrderedMap);

// Function to decompose a string with snarl information
// Given a string representing a path, return a vector of edges
std::vector<std::string> decompose_string(const std::string& s);

// Function to determine and extract an integer from the string
// TODO: I think this takes a string representing a path, the length of the string, and an index and returns the substring of the node on the path at i
std::pair<int, std::string> determine_str(const std::string& s, size_t length_s, size_t i);

// Function to decompose a list of snarl strings
// Given a vector of paths, return the same vector with the paths represented as vectors of edges
const std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst);

// Given a path through the snarl, a matrix of edges for each sample/haplotype, and the number of columns (samples/haplotypes),
// return the columns for which all edges (rows) in the path are set, i.e. the haplotypes with the given path.
std::vector<size_t> identify_correct_path(const std::vector<std::string>& decomposed_snarl,
                                        const Matrix& matrix,
                                        const size_t num_cols);

std::vector<std::vector<size_t>> transpose_matrix(const std::vector<std::vector<size_t>>& matrix);

void find_two_largest_indices(const std::vector<size_t>& vec, size_t& major_index_1, size_t& major_index_2);

#endif
