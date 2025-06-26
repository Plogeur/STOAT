#ifndef LIST_SNARL_PATHS
#define LIST_SNARL_PATHS

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <iostream>
#include <array>
#include <chrono>
#include <cassert>
#include <regex>
#include <stdexcept>
#include <utility>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include <bdsg/overlays/packed_path_position_overlay.hpp>
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/path_handle_graph.hpp>
#include "utils.hpp"

using namespace std;
using namespace bdsg;
using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

// A class representing a path as a vector of strings representing nodes
// TODO: This is only used in fill_pretty_paths()
class Path {
private:
    std::vector<std::string> nodes;
    std::vector<char> orients;

public:
    // Constructor
    Path();

    // Add a node with known orientation
    void addNode(const std::string& node, char orient);

    // Add a node handle and extract information using the string representation
    bool addNodeHandle(const net_handle_t& node_h, const SnarlDistanceIndex& stree);

    // Get the string representation of the path
    std::string print() const;

    // Flip the path orientation
    void flip();

    // Get the size of the path
    size_t size() const;

    // Count the number of reversed nodes
    size_t nreversed() const;
};

// Load the distance index and graph and return unique_ptrs to them
std::tuple<std::unique_ptr<bdsg::SnarlDistanceIndex>, 
           std::unique_ptr<bdsg::PackedGraph>, 
           handlegraph::net_handle_t, 
           std::unique_ptr<bdsg::PackedPositionOverlay>>
parse_graph_tree(const std::string& pg_file, const std::string& dist_file);

// Function to calculate the type of variant
vector<string> calcul_pos_type_variant(const vector<tuple<string, size_t, size_t, size_t, size_t, bool>>& list_length_paths);

// Function to find snarl ID
string find_snarl_id(SnarlDistanceIndex& stree, net_handle_t& snarl);

// Function to follow edges
// Following the netgraph edges from the last element of path, and make a new path for each different continuation of the path.
// If the new path(s) reach the end of the snarl, add the new path to finished_paths. Otherwise, add it to paths.
void follow_edges(SnarlDistanceIndex& stree,
    vector<vector<net_handle_t>>& finished_paths,
    const vector<net_handle_t>& path,
    vector<vector<net_handle_t>>& paths,
    PackedGraph& pg,
    const bool& cycle);

// Function to save snarls
// Returns a vector of <snarl net handle, reference path name, start offset on reference, end offset on reference, does the reference pass through the snarl> 
// for each snarl in the snarl tree 
// If the reference doesn't pass through the snarl, then the snarl's reference path and offsets will be the same as its lowest ancestor that has reference coordinates. 
vector<tuple<net_handle_t, string, size_t, size_t, bool>> save_snarls(
                            SnarlDistanceIndex& stree, 
                            net_handle_t& root,
                            PackedGraph& pg, 
                            unordered_set<string>& ref_paths,
                            PackedPositionOverlay& ppo);

// Function to fill pretty paths
// Given a vector of paths finished_paths (as vectors of net_handle_ts of children of a snarl), return 
// a vector of paths and their corresponding variant types  
tuple<vector<string>, vector<string>> fill_pretty_paths(
                            SnarlDistanceIndex& stree, 
                            PackedGraph& pg, 
                            vector<vector<net_handle_t>>& finished_paths);

// Function to loop over snarls and write output to output_file
// Output is a tsv of <chromosome, start pos, end pos, snarl, paths, variant type, reference>
// Returns a map from chromosome name to a vector of <snarl name, paths, start position, end position, variant type>
std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>> loop_over_snarls_write(
                            SnarlDistanceIndex& stree, 
                            vector<tuple<net_handle_t, string, size_t, size_t, bool>>& snarls, 
                            PackedGraph& pg, 
                            const string& output_file, 
                            const string& output_snarl_not_analyse, 
                            const size_t& children_treshold,
                            const size_t& path_length_threshold,
                            const size_t& cycle_threshold,
                            bool bool_return);

#endif
