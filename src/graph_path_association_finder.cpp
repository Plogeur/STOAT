#include "graph_path_association_finder.hpp"
#include "utils.hpp"
#include "writer.hpp"
#include "binary_table.hpp"

//#define DEBUG_ASSOCIATION_FINDER

namespace stoat_graph {

AssociationFinder::AssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                                     const bdsg::SnarlDistanceIndex& distance_index,
                                     std::shared_ptr<Partitioner> partitioner,
                                     const std::set<std::string>& samples_of_interest, 
                                     std::string reference_sample,
                                     std::string test_method,
                                     size_t allele_size_limit,
                                     std::ostream& out_associated,
                                     std::ostream& out_unassociated) :
    graph(graph), 
    distance_index(distance_index), 
    partitioner(std::move(partitioner)),
    samples_of_interest(samples_of_interest), 
    reference_sample(std::move(reference_sample)),
    test_method(std::move(test_method)),
    allele_size_limit(allele_size_limit),
    out_associated(out_associated),
    out_unassociated(out_unassociated)
    {}

void AssociationFinder::test_snarls() const {

    //TODO: Make this general
    // If the file output has a header, write it ?
    // Matis ans : why just do an if binary/quantitative ?
    stoat_vcf::write_binary_header(out_associated);

    std::vector<handlegraph::net_handle_t> chains;
    chains.reserve(graph.get_node_count()/100);
    handlegraph::net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, [&] (handlegraph::net_handle_t chain) {
        chains.emplace_back(chain);
        return true;
    });

    FisherKhi2 fisher_chi2_tester;
    while (!chains.empty()) {
        handlegraph::net_handle_t chain = chains.back();
        chains.pop_back();

        distance_index.for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {

            //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
            if (distance_index.is_snarl(snarl) && snarl_is_eligible(snarl) ) {
#ifdef DEBUG_ASSOCIATION_FINDER
                cerr << "Test snarl " << distance_index.net_handle_as_string(snarl) << endl;
#endif

                // Should we write this?
                bool write_output = false;

                // the strings we are going to output
                string group_paths, allele_number_str, min_row_index_str, 
                        numb_colum_str, inter_group_str, average_str,
                        fastfisher_p_value, chi2_p_value = "NA";
                string variant_type = "UNKNOWN_TYPE";

                // Each set represents a partition of samples that takes the same path through the snarl's netgraph
                std::vector<std::set<std::string>> sample_partitions = partitioner->partition_samples_in_snarl(graph, distance_index, snarl);

                if (test_method == "exact") {
#ifdef DEBUG_ASSOCIATION_FINDER
                        cerr << "\tTRUTH" << endl;
                        for (const std::string& sample : samples_of_interest) {
                            cerr << "\t\t" << sample << endl;
                        }
#endif

                    for (const std::set<std::string>& partition : sample_partitions) {
#ifdef DEBUG_ASSOCIATION_FINDER
                        cerr << "\tPARTITION" << endl;
                        for (const std::string& sample : partition) {
                            cerr << "\t\t" << sample << endl;
                        }
#endif
                        if (partition == samples_of_interest) {
#ifdef DEBUG_ASSOCIATION_FINDER
                            cerr << "\tFound exact match" << endl;
#endif
                            write_output = true;
                            break;
                        }
                    }

                } else {

                    // If we are using a real statistical test, then always write the output because the BH correction will need all the p-values
                    // TODO: This could do what pangwas was doing to keep track of only good p-values instead of writing everything
                    write_output = true;

                    // Fill in the genotypes. Each item in these vectors is an allele (path/sample partition)
                    std::vector<size_t> genotype_associated(sample_partitions.size(), 0);
                    std::vector<size_t> genotype_unassociated(sample_partitions.size(), 0);
                    for (size_t i = 0 ; i < sample_partitions.size() ; i++) {
                        const std::set<std::string> sample_set = sample_partitions[i];
                        for (const std::string sample : sample_set) {
                            if (samples_of_interest.count(sample) != 0) {
                                genotype_associated[i]++;
                            } else {
                                genotype_unassociated[i]++;
                            }
                        }
                    }

                    //Get a bunch of strings that get used for the output
                    // TODO: This function should probably be part of the output function
                    auto [group_paths, 
                        allele_number_str, min_row_index_str, 
                        numb_colum_str, inter_group_str, average_str] = stoat_vcf::binary_stat_test(genotype_associated, genotype_unassociated);
 
                    // Run the statistical test
                    auto [fastfisher_p_value, chi2_p_value] = fisher_chi2_tester.fisher_khi2(genotype_associated, genotype_unassociated);

                }
                // TODO idk what to put for chr
                string chr = "NA"; 
                // Matis ans : why don't you put the actual chr ref if the snarl containt it and something like not_ref if it's not
                //TODO: Maybe I sould keep the snarls as snarl_data_t's? 
                // TODO: get the type properly
                stoat_vcf::Snarl_data_t snarl_data_s(snarl, graph, distance_index);

                // Get the offsets of the start and end nodes along the reference
                std::vector<path_range_t> ranges = get_coordinates_of_snarl(graph, distance_index, snarl, true, "", false);
                if (ranges.size() != 0) {
                    snarl_data_s.start_positions = graph.get_position_of_step(ranges.front().start);
                    snarl_data_s.end_positions = graph.get_position_of_step(ranges.front().end);

                    chr = graph.get_path_name(graph.get_path_handle_of_step(ranges.front().start));
                }
                
                if (write_output) {
                    # pragma omp critical (out_associated) 
                    {
                        stoat_vcf::write_binary(out_associated, chr, snarl_data_s, variant_type, fastfisher_p_value, chi2_p_value, "", allele_number_str, min_row_index_str,
                                     numb_colum_str, inter_group_str, average_str, group_paths);
                    }
                }

                // Add the child chains to the stack
                distance_index.for_each_child(snarl, [&] (handlegraph::net_handle_t child) {

                    chains.emplace_back(child);
                    return true;
                });
            }
            return true;
        });
    }
}

bool AssociationFinder::snarl_is_eligible(const handlegraph::net_handle_t& snarl) const {
    return distance_index.maximum_length(snarl) >= allele_size_limit;
}

}//end pangwas namespace

