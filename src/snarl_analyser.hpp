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

class SnarlAnalyser {
    public:
        SnarlAnalyser();
        virtual ~SnarlAnalyser()=default;

};

class BinaryAnalyser : public SnarlAnalyser {
    public:
        BinaryAnalyser();
        ~BinaryAnalyser()=default;

};

class QuantitativeAnalyser : public SnarlAnalyser {
    public:
        QuantitativeAnalyser();
        ~QuantitativeAnalyser()=default;

};

class EqtlAnalyser : public SnarlAnalyser {
    public:
        EqtlAnalyser();
        ~EqtlAnalyser()=default;

};


/// Given a snarl_data_s for one snarl, make a genotype matrix and write the tsv output
void write_snarl_line_binary(const EdgeBySampleMatrix& edge_matrix, 
    const Snarl_data_t& snarl_data_s,
    const std::vector<bool>& binary_phenotype, 
    const std::string& chr,
    const std::vector<std::vector<double>>& covar,
    const double& maf,  
    const double& table_threshold, 
    const std::string& output_dir, 
    size_t sample_count,
    std::ofstream& outf);

/// Similar to write_snarl_line_binary, get the genotypes and write the tsv output
void write_snarl_line_quantitative(const EdgeBySampleMatrix& edge_matrix, 
    const Snarl_data_t& snarl_data_s,
    const std::vector<double>& quantitative_phenotype, 
    const std::string &chr,
    const std::vector<std::vector<double>>& covar,
    const double& maf,  
    const double& table_threshold, 
    const std::string& output_dir, 
    size_t sample_count,
    std::ofstream& outf);

/// Similar to write_snarl_line_binary and write_snarl_line_quantitative, get the genotype and write the tsv output
void write_snarl_line_eqtl(const EdgeBySampleMatrix& edge_matrix, 
    const Snarl_data_t& snarl_data_s,
    const std::vector<stoat_vcf::Qtl_data>& eqtl,
    const std::string& chr, 
    const std::vector<std::vector<double>>& covar,
    const double& maf,  
    const double& table_threshold, 
    const std::string& regression_dir,
    const size_t& windows_gene_threshold, 
    size_t sample_count,
    std::ofstream& outf);

/// Given a list of paths in a snarl, return vectors of counts of each sample taking an allele for the two alleles with the highest counts over all samples
/// Each vector returned corresponds to one allele, each entry in the vector is a sample, the value is a count of the number of times a sample takes the path/allele (probably binary?)
std::pair<std::vector<size_t>, std::vector<size_t>> create_table_short_path(const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl, size_t sample_count, const EdgeBySampleMatrix& edge_matrix);

/// For each snarl, write a bim file and a bed file (PLINK formats)
void create_bim_bed(const std::vector<Snarl_data_t>& snarls, size_t sample_count, 
    const EdgeBySampleMatrix& edge_matrix,
    std::string chromosome, std::ofstream& outbim, std::ofstream& outbed);

/// Return true if any column exceeds the MAF threshold 
bool check_MAF_threshold_quantitative(const std::vector<std::vector<double>>& df, const double& maf);

bool check_MAF_threshold_binary(
    const std::vector<size_t>& g0, const std::vector<size_t>& g1,
    const size_t& totalSum, const size_t& length_column_headers, 
    const double& maf);

/// Go through the vcf by chromosome, parse it to get a matrix of genotypes (either binary, quantitative, or eqtl, depending on the phenotype type),
/// then write the output (also depending on the phenotype type).
/// window_gene_threshold and eqtl_map are only used for eqtl output 
void chunk_chromosome_and_write_tsv(phenotype_type_t phenotype_type,
     htsFile* &ptr_vcf,
     bcf_hdr_t* &hdr,
     bcf1_t* &rec,
     const std::vector<std::string> &list_samples,
     const std::unordered_map<std::string, std::vector<Snarl_data_t>> &chr_to_snarl_data,
     const std::vector<bool>& binary_pheno,
     const std::vector<double>& quantitative_pheno,
     const std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>>& eqtl_map,
     std::vector<std::vector<double>> covar,
     const double& maf,
     const double& table_threshold,
     const size_t& windows_gene_threshold,
     const std::string& regression_dir,
     const std::string& output_filename);

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& snarl_chr,
    const std::string& output_dir);

std::vector<size_t> found_gene_snarl(
    const std::vector<Qtl_data>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold);

void create_fam(const std::vector<std::pair<std::string, int>> &pheno, 
    const std::string& output_path);

/// Make an EdgeBySampleMatrix representing the genotypes in a vcf and the pointers to the vcf but advanced to the end of the chromosome?
std::tuple<EdgeBySampleMatrix, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const std::vector<std::string>& sample_names, std::string &chr, size_t &num_paths_ch);

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
