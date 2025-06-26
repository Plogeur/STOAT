#ifndef snarl_data_t
#define snarl_data_t

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

struct Snarl_data_t {

    public:
        // Maximum number of snarls to reserve space for
        void reserve(size_t snarl_count = MAX_SNARLS);

        // Add a snarl
        void add_snarl(const std::pair<size_t, size_t>& name, const Path_traversal_t& paths,
                    size_t start, size_t end, const std::vector<std::string>& path_nodes);

    private:

        std::vector<std::pair<size_t, size_t>> snarl_id;
        std::vector<Path_traversal_t> snarl_paths;
        std::vector<size_t> start_positions;
        std::vector<size_t> end_positions;
        std::vector<std::vector<std::string>> type_variants; // because of complexe X/Y we can't use size_t here
};

struct Node_traversal_t {
    private:
        size_t node_id : 63;
        bool is_reverse : 1;

    public:
        Node_traversal_t(const size_t &id, const bool &rev);
        
        // Convert to string representation
        std::string to_string() const;
};

struct Path_traversal_t {
    private:
        std::vector<Node_traversal_t> paths; // Nodes in the path

    public:
        // add a node traversal to the path
        Path_traversal_t() = default;
        add_node_traversal_t(const Node_traversal_t &paths);

        // convert to string representation
        std::string to_string() const;
};

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
    Path_traversal_t print() const;

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
vector<string> calcul_pos_type_variant(const vector<tuple<string, size_t, size_t, size_t, size_t, bool>>& list_length_paths);

// Function to find snarl ID
std::pair<size_t, size_t> find_snarl_id(SnarlDistanceIndex& stree, net_handle_t& snarl);

// Function to follow edges
void follow_edges(SnarlDistanceIndex& stree,
    vector<vector<net_handle_t>>& finished_paths,
    const vector<net_handle_t>& path,
    vector<vector<net_handle_t>>& paths,
    PackedGraph& pg,
    const bool& cycle);

// Function to save snarls
vector<tuple<net_handle_t, string, size_t, size_t, bool>> save_snarls(
                            SnarlDistanceIndex& stree, 
                            net_handle_t& root,
                            PackedGraph& pg, 
                            unordered_set<string>& ref_paths,
                            PackedPositionOverlay& ppo);

// Function to fill pretty paths
tuple<vector<Path_traversal_t>, vector<string>> fill_pretty_paths(
                            SnarlDistanceIndex& stree, 
                            PackedGraph& pg, 
                            vector<vector<net_handle_t>>& finished_paths);

// Function to loop over snarls and write output
std::unordered_map<std::string, Snarl_data_t> loop_over_snarls_write(
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
