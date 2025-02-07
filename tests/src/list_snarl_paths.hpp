#ifndef LIST_SNARL_PATHS
#define LIST_SNARL_PATHS

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_map>
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

using namespace std;
using namespace bdsg;
using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

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
    void addNodeHandle(const net_handle_t& node_h, const SnarlDistanceIndex& stree);

    // Get the string representation of the path
    std::string print() const;

    // Flip the path orientation
    void flip();

    // Get the size of the path
    size_t size() const;

    // Count the number of reversed nodes
    size_t nreversed() const;
};

std::tuple<std::unique_ptr<bdsg::SnarlDistanceIndex>, 
           std::unique_ptr<bdsg::PackedGraph>, 
           handlegraph::net_handle_t, 
           std::unique_ptr<bdsg::PackedPositionOverlay>>
parse_graph_tree(const std::string& pg_file, const std::string& dist_file);

// Function to calculate the type of variant
pair<vector<string>, size_t> calcul_pos_type_variant(const vector<vector<int>>& list_list_length_paths);

// Function to check threshold
void check_threshold(double proportion);

// Function to find snarl ID
string find_snarl_id(SnarlDistanceIndex& stree, net_handle_t& snarl);

// Function to follow edges
void follow_edges(
    SnarlDistanceIndex& stree,
    vector<vector<net_handle_t>>& finished_paths,
    vector<net_handle_t>& path,
    vector<vector<net_handle_t>>& paths,
    PackedGraph& pg
);

// Function to save snarls
vector<tuple<net_handle_t, string, size_t>> save_snarls(
                            SnarlDistanceIndex& stree, 
                            net_handle_t& root,
                            PackedGraph& pg, 
                            unordered_set<string>& ref_paths,
                            PackedPositionOverlay& ppo);

// Function to fill pretty paths
tuple<vector<string>, vector<string>, size_t> fill_pretty_paths(
                            SnarlDistanceIndex& stree, 
                            PackedGraph& pg, 
                            vector<vector<net_handle_t>>& finished_paths);

// Function to loop over snarls and write output
pair<vector<tuple<string, vector<string>, string, string, string>>, int> loop_over_snarls_write(
                            SnarlDistanceIndex& stree, 
                            vector<tuple<net_handle_t, string, size_t>>& snarls, 
                            PackedGraph& pg, 
                            const string& output_file, 
                            const string& output_snarl_not_analyse, 
                            int children_treshold, 
                            bool bool_return);

#endif 
