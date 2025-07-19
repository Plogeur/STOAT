#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"
#include "snarl_data_t.hpp"

using namespace std;
namespace stoat_vcf{

// Write headers
void write_binary_header(std::ofstream& outstream);
void write_binary_covar_header(std::ofstream& outstream);
void write_quantitative_header(std::ofstream& outstream);
void write_eqtl_header(std::ofstream& outstream);

// Write lines
void write_binary_covar(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                        const std::string& beta, const std::string& se, size_t allele_number, const std::vector<size_t>& allele_paths);

void write_eqtl(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                    const std::string& gene_name, const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                    const std::string& beta, const std::string& se, size_t allele_number, const std::vector<size_t>& allele_paths);

void write_binary(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                    const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& p_value_adjusted, 
                    const std::string& allele_number_str, const std::string& min_row_index_str, const std::string& num_colum_str,
                    const std::string& inter_group_str, const std::string& average_str, const std::string& group_paths);

void write_quantitative(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                        const std::string& beta, const std::string& se, size_t allele_number, const std::vector<size_t>& allele_paths);

void write_fasta(std::ofstream& outstream, const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                        const handlegraph::net_handle_t& snarl, const std::unordered_map<std::string, bool>& samples, const string& reference_name);

void writeSignificantTableToTSV(
    const std::vector<std::vector<double>>& table,
    const std::vector<std::string>& list_snarl,
    const std::vector<std::string>& list_samples,
    const std::string& filename);

} //end namespace

#endif
