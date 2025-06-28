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
#include <cstddef>
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

// Container for snarl data : 
// const std::pair<size_t, size_t>& snarl_id_,
// const std::vector<Path_traversal_t>& snarl_paths_,
// size_t start_positions_, size_t end_positions_,
// const std::vector<std::string>& type_variants_

struct Node_traversal_t { // 64 bits per node 
    private:
        size_t node_id : 63; // 63 bits for node ID
        bool is_reverse : 1; // 1 bit for orientation (true for reverse, false for forward)

    public:
        Node_traversal_t(const size_t &id, const bool &rev);
        
        // Getters
        size_t get_node_id() const;
        bool get_is_reverse() const;

        // Convert to string representation
        std::string to_string() const;

        bool operator==(const Node_traversal_t& other) const;
};

struct Edge_t { // 128 bits per edge 
    private:
        std::pair<Node_traversal_t, Node_traversal_t> edge;

    public:
        Edge_t(const Node_traversal_t &node_traversal_1, const Node_traversal_t &node_traversal_2);
        
        // Converter
        std::pair<size_t, size_t> print_pair_node() const;

        // Accessor to edge, useful for hashing and comparison
        const std::pair<Node_traversal_t, Node_traversal_t>& get_edge() const;

        // Comparison operator
        bool operator==(const Edge_t &other) const;
};

namespace std {
    template <>
    struct hash<Node_traversal_t> {
        size_t operator()(const Node_traversal_t& node) const {
            // Simple way: Shift node_id and pack is_reverse into the lower bit
            return (node.get_node_id() << 1) | static_cast<size_t>(node.get_is_reverse());
        }
    };
}

namespace std {
    template <>
    struct hash<Edge_t> {
        size_t operator()(const Edge_t& edge) const {
            const auto& pair = edge.get_edge();
            size_t h1 = hash<Node_traversal_t>()(pair.first);
            size_t h2 = hash<Node_traversal_t>()(pair.second);
            
            // Standard hash combination
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}

struct Path_traversal_t {
    private:
        std::vector<Node_traversal_t> paths; // Nodes in the path

    public:
        // add a node traversal to the path
        Path_traversal_t() = default;
        void add_node_traversal_t(const Node_traversal_t &paths);

        // Getters
        const std::vector<Node_traversal_t>& get_paths() const;
        
        // convert to string representation
        std::string to_string() const;
};

struct Snarl_data_t {
    public:
        // Constructor definition
        Snarl_data_t(const std::pair<size_t, size_t>& snarl_id_,
                    const std::vector<Path_traversal_t>& snarl_paths_,
                    const size_t start_positions_, const size_t end_positions_,
                    const std::vector<std::string>& type_variants_);  // Assuming path_nodes correspond to type_variants

        // Getters
        const std::pair<size_t, size_t>& get_snarl_id() const;
        const std::vector<Path_traversal_t>& get_snarl_paths() const;
        const size_t& get_start_positions() const;
        const size_t& get_end_positions() const;
        const std::vector<std::string>& get_type_variants() const;
        const std::tuple<std::string, std::vector<Path_traversal_t>, size_t, size_t, std::vector<std::string>>& get_snarl() const;

    private:
        std::vector<std::string> type_variants;
        std::vector<Path_traversal_t> snarl_paths;
        std::pair<size_t, size_t> snarl_id;
        size_t start_positions;
        size_t end_positions;
};

// Converter
std::string pairToString(const std::pair<size_t, size_t>& name);
std::pair<size_t, size_t> stringToPair(const std::string& str);
std::string vectorPathToString(const std::vector<Path_traversal_t>& vec_paths);
std::vector<Path_traversal_t> stringToVectorPath(std::string& str);

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
std::unordered_map<std::string, std::vector<Snarl_data_t>> loop_over_snarls_write(
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
