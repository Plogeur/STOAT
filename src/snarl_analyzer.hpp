#ifndef snarl_analyzer_HPP
#define snarl_analyzer_HPP

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
#include "stats_test.hpp"

using namespace std;

namespace stoat_vcf {

class SnarlAnalyzer {
public:
    SnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
                  const std::vector<std::vector<double>>& covariate, double maf, double table_threshold);

    ~SnarlAnalyzer()=default;

    /// Go through the vcf by chromosome, parse it to get a matrix of genotypes (either binary, quantitative, or eqtl, depending on the phenotype type),
    /// then write the output (also depending on the phenotype type).
    /// This calls write_header() to write the appropriate output header and analyze_and_write_snarl() for each snarl
    void process_snarls_by_chromosome_chunk(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec,
                                            const std::string& regression_dir, const std::string& output_filename);

    /// Update the EdgeBySampleMatrix representing the genotypes in a vcf and the pointers to the vcf but advanced to the end of the chromosome?
    std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> make_edge_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, std::string &chr, size_t &num_paths_ch);

    /// For the given snarl, analyze the snarl and write it to outf
    virtual void analyze_and_write_snarl(const std::string& chr, const Snarl_data_t& snarl_data, const std::string& regression_dir, std::ofstream& outf) = 0;

    /// Write the header of the output tsv file
    /// This should ideally call a write_header() function from writer.hpp to keep things consistent
    virtual void write_header(std::ofstream&outf) = 0;

//////////////// Private data members
protected:
    
    // Map chromosome name to a vector of snarl_data_t
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data;

    // A list of sample names
    const std::vector<std::string>& list_samples;

    // Covariate matrix
    const std::vector<std::vector<double>>& covariate;

    // Matrix of edges in each sample/haplotype
    // This generally is a per-chromosome or per-chunk matrix, so it must be updated for each new chunk being analyzed 
    EdgeBySampleMatrix& edge_matrix;
    const double& maf; 
    const double& table_threshold;
    std::ofstream& outf;
    std::string& chr; 
    const std::string& regression_dir;
    const std::unordered_map<std::vector<stoat_vcf::Qtl_data>>& eqtl_map;
    const size_t& windows_gene_threshold;
};

class BinarySnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinarySnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
                  const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, const std::vector<bool>& binary_phenotype);

    void analyze_and_write_snarl(const std::string& chr, const Snarl_data_t& snarl_data, const std::string& regression_dir, std::ofstream& outf);

    void write_header(std::ofstream&outf);

/////////////////// Private data members
protected:

    const std::vector<bool>& binary_phenotype;
    FisherKhi2& fk;
};

class BinaryCovarSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    BinaryCovarSnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
                  const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, const std::vector<bool>& binary_phenotype);

    void analyze_and_write_snarl(const std::string& chr, const Snarl_data_t& snarl_data, const std::string& regression_dir, std::ofstream& outf);

    void write_header(std::ofstream&outf);

/////////////////// Private data members
protected:

    const std::vector<bool>& binary_phenotype;
    LogisticRegression& lr;
};

class QuantitativeSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    QuantitativeSnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
                  const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, const std::vector<double>& quantitative_phenotype);

    void analyze_and_write_snarl(const std::string& chr, const Snarl_data_t& snarl_data, const std::string& regression_dir, std::ofstream& outf) ;

    void write_header(std::ofstream&outf);

/////////////////// Private data members
protected:

    const std::vector<double>& quantitative_phenotype;
    LinearRegression& lr;
};

class EQTLSnarlAnalyzer : public SnarlAnalyzer {

public:
    
    EQTLSnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
                  const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, 
                  const std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>>& eqtl_map,
                  size_t windows_gene_threshold);

    void analyze_and_write_snarl(const std::string& chr, const Snarl_data_t& snarl_data, const std::string& regression_dir, std::ofstream& outf);

    void write_header(std::ofstream&outf);

/////////////////// Private data members
protected:

    // TODO idk what these are 
    // Maps something to something else?
    // Matis ans : eqtl_map is an {chr name : std::vector<stoat_vcf::Qtl_data>}
    // is organise like that in the first place to optimize edge_matrix / eqtl linking
    // but now we can just use std::vector<stoat_vcf::Qtl_data> because we already know the chr
    // that we gonna use
    const std::vector<stoat_vcf::Qtl_data>& eqtl;
    const size_t& windows_gene_threshold;
    LinearRegression& lr;
};

/// Given a list of paths in a snarl, return vectors of counts of each sample taking an allele for the two alleles with the highest counts over all samples
/// Each vector returned corresponds to one allele, each entry in the vector is a sample, the value is a count of the number of times a sample takes the path/allele (probably binary?)
std::pair<std::vector<size_t>, std::vector<size_t>> create_table_short_path(
    const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl, 
    size_t sample_count, 
    const EdgeBySampleMatrix& edge_matrix);

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
