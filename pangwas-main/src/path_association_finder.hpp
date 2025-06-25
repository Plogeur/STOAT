#ifndef PANGWAS_PATH_ASSOCIATION_FINDER_HPP_INCLUDED
#define PANGWAS_PATH_ASSOCIATION_FINDER_HPP_INCLUDED

#include "association_finder.hpp"
#include "utils.hpp"

using namespace std;
namespace pangwas{

/*
    A class for finding associations based on the paths in the graph
*/
class PathAssociationFinder : public AssociationFinder {
     public:

        /// Initialize as with the base class, and also find the set of all paths in the graph
        PathAssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                              const bdsg::SnarlDistanceIndex& distance_index,
                              std::string test_meethod,
                              const std::set<std::string>& samples_of_interest,
                              std::string reference_name, 
                              std::string output_format,
                              std::ostream& out_associated,
                              std::ostream& out_unassociated,
                              size_t allele_size_limit,
                              double p_value);


     protected:

        /// Get partitions of samples in the snarl 
        std::vector<std::set<std::string>> partition_samples_in_snarl(const handlegraph::net_handle_t& snarl) const;

        /// Given a snarl, partition the paths going through the snarl based on the walks they take in the netgraph.
        /// Unlike get_start_edge_sets, any path not in the snarl will also be returned as a separate set
        /// Returns sets of samples + haplotypes
        /// TODO: Maybe it shouldn't but I included it since there may be tips
        std::vector<std::set<sample_hap_t>> get_walk_sets(const bdsg::net_handle_t& snarl) const;

        /// Given a snarl, partition the paths going through the snarl based on the edges going into the snarl from the start bound.
        /// If a path traverses the snarl multiple times, it may appear in multiple sets
        /// Returns sets of sample + haplotypes
        /// TODO: I'm also not sure if this is the correct behavior
        std::vector<std::set<sample_hap_t>> get_start_edge_sets(const bdsg::net_handle_t& snarl) const;

    protected: 
        // Additional members

        /// A set of all samples+haplotypes in the graph
        std::set<sample_hap_t> all_sample_haplotypes;

};

}

#endif
