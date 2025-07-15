#ifndef PANGWAS_ASSOCIATION_FINDER_HPP_INCLUDED
#define PANGWAS_ASSOCIATION_FINDER_HPP_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "tester.hpp"

using namespace std;
namespace pangwas{

/***
    General template class for finding associations in a graph.
    This will implement the helper functions needed to traverse the graph, filter snarls that 
    are too small, write the output, etc.
    Inherited classes must implement is_snarl_associated() to test each snarl.
***/
class AssociationFinder {

    protected:
        // At a minimum, an AssociationFinder must have a graph (with path information for printing reference coordinates),
        // a distance index, the names of the samples we are interested in, and, optionally, the name of the reference
        const handlegraph::PathPositionHandleGraph& graph;
        const bdsg::SnarlDistanceIndex& distance_index; 
        const std::string reference_name;
        const std::set<std::string>& samples_of_interest;
        const std::string output_format;
        std::ostream& out_associated = std::cout;
        std::ostream& out_unassociated = std::cout;
        size_t allele_size_limit;

        // object for doing the actual association test
        std::shared_ptr<Tester> tester;


    public:

        /// Create an association finder with the graph and distance index, a set of samples of interest for which we want associated
        /// variants, a std::string of the reference sample name (may be empty), the output format (tsv or fasta), 
        /// filenames for writing associated alleles and unassociated alleles, and a size limit for the minimum length of snarl reported,
        /// measured as the "maximum" length of a snarl
        /// 
        /// This doesn't fill in the tester. It is the derived class's responsibility because otherwise it would duplicate a lot of work
        AssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                          const bdsg::SnarlDistanceIndex& distance_index, 
                          std::string test_method, 
                          const std::set<std::string>& samples_of_interest,  
                          std::string reference_name,
                          std::string output_format, std::ostream& out_associated, 
                          std::ostream& out_unassociated,
                          size_t allele_size_limit, double p_value);

        
        /// Main function that gets called to go through the graph, call is_snarl_associated, and
        /// write the output. 
        void write_associated_snarls() const;

    protected:
        //////////////////////////////////  Template functions for determining association

        /// Given a snarl, return a partition of samples that will be used to determine association with samples of interest.
        /// This is a template function that must be implemented by inherited classes
        virtual std::vector<std::set<std::string>> partition_samples_in_snarl(const handlegraph::net_handle_t& snarl) const = 0;
        

    protected:
        /////////////////////////// Non-template functions for determining association 

        /// Is the given snarl associated with the trait shared by samples_of_interest?
        /// If yes, return true and also a set of sample names to be output (one per representative partition for unassociated samples)
        /// This calls partition_samples_in_snarl then uses the Tester to test the significance
        std::pair<bool, std::unordered_set<std::string>> is_snarl_associated(const handlegraph::net_handle_t& snarl) const;


    protected:

        
        ////////////////////////////////// Functions for writing output

        // Do we care about this snarl? Based on allele_size_limit
        bool snarl_is_eligible(const handlegraph::net_handle_t& snarl) const;

        // Write the header of a file, depending on output_format 
        void write_header() const;

        // Write a single snarl
        void write_snarl(const handlegraph::net_handle_t& snarl, const std::unordered_set<std::string>& samples) const;

        // Given a snarl, output the coordinates of the snarl as a tsv of
        // reference path, start offset, end offset, max length
        void write_tsv_of_snarl(const handlegraph::net_handle_t& snarl) const; 

        // Given a snarl, output a fasta of all paths going through the snarl
        // If sample names are given, only output the records for the given samples 
        void write_fasta_of_snarl(const handlegraph::net_handle_t& snarl, const std::unordered_set<std::string>& samples) const;

        ///////////////////////////////// Helper functions doing stuff on the graph

        // A struct for holding a range along the path
        struct path_range_t {
            handlegraph::step_handle_t start;
            handlegraph::step_handle_t end;
        };

        /// Given a snarl, return a vector of path_ranges of that snarl (the boundary nodes).
        /// Since a path can traverse a snarl multiple times, this returns each start-to-end (or end-to-start) range
        /// of step_handle_t's, ordered according to the order of the path.
        /// If the path leaves by the same bound (for example start-> start<- start-> end->), then the range will include
        /// the outermost start->end range.
        /// If get_reference is true, return a reference path and its coordinates.
        /// This will first try to find a path with the sample name, if not empty, then a reference-sense path, then with any path traversing the snarl.
        /// If get_reference is false, try to find coordinates on a path containing the given sample name, or if it fails, with any path. 
        /// If get_reference is false and sample_name is empty and get_all_paths is true, return all coordinates for all paths
        std::vector<stoat::path_range_t> get_coordinates_of_snarl(const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths) const;

        /// The function that gets called by get_coordinates_of_snarl
        /// This either looks for a particular sample, or a reference-sense path, or all paths
        std::vector<stoat::path_range_t> stoat::get_coordinates_of_snarl(const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths) const;

};


}

#endif
