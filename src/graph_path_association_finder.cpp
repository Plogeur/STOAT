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
                                     std::string test_method,
                                     size_t total_sample_count,
                                     size_t allele_size_limit,
                                     std::ofstream& out_associated,
                                     std::ofstream& out_unassociated) :
    graph(graph), 
    distance_index(distance_index), 
    partitioner(std::move(partitioner)),
    samples_of_interest(samples_of_interest), 
    test_method(std::move(test_method)),
    total_sample_count(total_sample_count),
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

    stoat_vcf::FisherKhi2 fk();
    while (!chains.empty()) {
        handlegraph::net_handle_t chain = chains.back();
        chains.pop_back();

        distance_index.for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {

            //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
            if (distance_index.is_snarl(snarl) && snarl_is_eligible(snarl) ) {

                std::vector<std::set<std::string>> sample_partitions = partitioner->partition_samples_in_snarl(graph, distance_index, snarl);

                if (test_method == "exact") {
                    // TODO add exact test here I supposed

                } else {
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

                    const auto& [group_paths, 
                        allele_number_str, min_row_index_str, 
                        numb_colum_str, inter_group_str, average_str] = 
                        stoat_vcf::binary_stat_test(genotype_associated, genotype_unassociated, 
                                group_paths, allele_number_str, min_row_index_str,
                                    numb_colum_str, inter_group_str, average_str);
                    
                    const auto& [fastfisher_p_value, chi2_p_value] = fk.fisher_khi2(g0, g1);

                    # pragma omp critical (out_associated) 
                    {
                        // TODO idk what to put for chr
                        // Matis ans : why don't you put the actual chr ref if the snarl containt it and something like not_ref if it's not
                        stoat_vcf::write_binary(out_associated, "?", snarl_data_s, type_var_str, fastfisher_p_value, chi2_p_value, "", allele_number_str, min_row_index_str,
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

