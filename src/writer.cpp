#include "writer.hpp"

//#define DEBUG_WRITER

namespace stoat_vcf {

void write_binary_covar_header(std::ofstream& outstream) {
    outstream << "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS" << endl;
}

void write_binary_header(std::ofstream& outstream) {
    outstream << "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tP_ADJUSTED\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS" << endl;
}

void write_quantitative_header(std::ofstream& outstream) {
    outstream << "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS" << endl;
}

void write_eqtl_header(std::ofstream& outstream) {
    outstream <<  "CHR\tPOS\tSNARL\tTYPE\tGENE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS" << endl;
}

void write_eqtl(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                   const std::string& gene_name, const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                   const std::string& beta, const std::string& se, size_t allele_number, const std::vector<size_t>& allele_paths) {
    outstream << chr << "\t" 
              << snarl_data_s.start_positions << "\t" 
              << pairToString(snarl_data_s.snarl_ids) << "\t" 
              << type_var_str << "\t" 
              << gene_name << "\t" 
              << p_value  << "\t" 
              << p_value_adjusted << "\t" 
              << r2 << "\t" 
              << beta << "\t" 
              << se << "\t" 
              << allele_number << "\t" 
              << vectorToString(allele_paths) << endl;

}

void write_binary_covar(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                        const std::string& beta, const std::string& se, size_t allele_number, const std::vector<size_t>& allele_paths) {
    outstream << chr << "\t" 
              << snarl_data_s.start_positions << "\t" 
              << pairToString(snarl_data_s.snarl_ids) << "\t" 
              << type_var_str << "\t" 
              << p_value << "\t" 
              << p_value_adjusted << "\t" 
              << r2 << "\t" 
              << beta << "\t" 
              << se << "\t" 
              << allele_number << "\t" 
              << vectorToString(allele_paths) << endl;

}

void write_binary(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& fastfisher_p_value, const std::string& chi2_p_value, const std::string& p_value_adjusted, 
                        const std::string& allele_number_str, const std::string& min_row_index_str, const std::string& num_colum_str,
                        const std::string& inter_group_str, const std::string& average_str, const std::string& group_paths) {
    outstream << chr << "\t" 
              << snarl_data_s.start_positions << "\t" 
              << pairToString(snarl_data_s.snarl_ids) << "\t" 
              << type_var_str << "\t" 
              << fastfisher_p_value << "\t" 
              << chi2_p_value << "\t" 
              << p_value_adjusted << "\t" 
              << allele_number_str << "\t" 
              << min_row_index_str << "\t" 
              << num_colum_str << "\t" 
              << inter_group_str << "\t" 
              << average_str << "\t" 
              << group_paths << endl;
}


void write_quantitative(std::ofstream& outstream, const std::string& chr, const Snarl_data_t& snarl_data_s, const std::string& type_var_str,
                        const std::string& p_value, const std::string& p_value_adjusted, const std::string& r2,
                        const std::string& beta, const std::string& se, size_t allele_number, const std::vector<size_t>& allele_paths) {
    outstream << chr << "\t" 
              << snarl_data_s.start_positions << "\t" 
              << pairToString(snarl_data_s.snarl_ids) << "\t" 
              << type_var_str << "\t" 
              << p_value  << "\t" 
              << p_value_adjusted << "\t" 
              << r2 << "\t" 
              << beta << "\t" 
              << se << "\t" 
              << allele_number << "\t" 
              << vectorToString(allele_paths) << "\n";

}
}//end namespace

