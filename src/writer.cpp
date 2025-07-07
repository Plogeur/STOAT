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
void write_fasta(std::ofstream& outstream, const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                        const handlegraph::net_handle_t& snarl, const std::unordered_map<std::string, bool>& samples, const string& reference_name) {
    
    // A handle_t of the start bound facing in
    handlegraph::handle_t start_handle = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph);
    
    // Get a unique name for the snarl, as the start and end ids
    // I think even in the case of a looping chain where there is another snarl on the other end, the order of node ids will be flipped
    std::string snarl_name = "snarl:" +
                             std::to_string((int)graph.get_id(start_handle)) +
                             "-" +
                             std::to_string((int)distance_index.node_id(distance_index.get_bound(snarl, true, false)));
    
    // Get a reference range for the snarl.
    // If the reference goes through the snarl multiple times, get the largest interval
    
    std::vector<path_range_t> ref_ranges = get_coordinates_of_snarl(graph, distance_index, snarl, true, reference_name, false);
    std::string ref_coordinates  = "NOREF:?:?";
    int start_offset = std::numeric_limits<int>::max();
    int end_offset = 0;
    //Only get the coordinates for one path
    bool first = true;
    handlegraph::path_handle_t ref_path;
    for (const path_range_t& ref_range : ref_ranges){
        if (first) {
            first = false;
            ref_path = graph.get_path_handle_of_step(ref_range.start);
        } else if (graph.get_path_handle_of_step(ref_range.start) != ref_path) {
            continue;
        }
        ref_coordinates = graph.get_path_name(graph.get_path_handle_of_step(ref_range.start));
        start_offset = std::min(start_offset,
                                (int)(graph.get_position_of_step(ref_range.start) +
                                     distance_index.minimum_length(distance_index.get_net(graph.get_handle_of_step(ref_range.start), &graph))));
        end_offset = std::max(end_offset, (int)graph.get_position_of_step(ref_range.end));
    }
    
    if (ref_ranges.size() != 0) {
        ref_coordinates += ":" + std::to_string(start_offset) + "-" + std::to_string(end_offset);
    }
    
    // Now go through each path that goes through the snarl and print the sequence
    std::vector<path_range_t> path_ranges = get_coordinates_of_snarl(graph, distance_index, snarl, false, "", true);
    for (const path_range_t& path_range : path_ranges) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(path_range.start);
        if (samples.empty() || samples.count(get_sample_name_from_path(graph, path)) != 0) {
            //If we aren't checking samples, or if this is a sample we want
    
            // Print the header
            outstream << ">" << snarl_name << "|"
                << ref_coordinates << "|"
                << graph.get_path_name(path) << ":"
                << (graph.get_position_of_step(path_range.start) + distance_index.minimum_length(distance_index.get_net(graph.get_handle_of_step(path_range.start), &graph))) << "-"    
                << graph.get_position_of_step(path_range.end) << endl;
    
            // Now print the sequence in 80bp chunks.
            // Keep a buffer to print 80 bp at a time
            std::string sequence_buffer = "";
            handlegraph::step_handle_t next_step = graph.get_next_step(path_range.start);
            while (next_step != path_range.end) {
                std::string node_seq = graph.get_sequence(graph.get_handle_of_step(next_step));
                while (node_seq.size() != 0) {
    
                    // Fill in sequence_buffer up to 80 characters
                    size_t to_add = 80 - sequence_buffer.size();
                    sequence_buffer += node_seq.substr(0, to_add);
                    node_seq.erase(0, to_add);
    
                    // If the buffer is full, write it and clear it
                    if (sequence_buffer.size() == 80) {
                        outstream << sequence_buffer << endl;
                        sequence_buffer.clear();
                    }
                }
                handlegraph::step_handle_t step = next_step;
                if (!graph.has_next_step(step)) {
                    break;
                }
                next_step = graph.get_next_step(step);
            }
            if (!sequence_buffer.empty()) {
                outstream << sequence_buffer << endl;
            }
        }
    }
}

}//end namespace

