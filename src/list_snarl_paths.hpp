#ifndef LIST_SNARL_PATHS
#define LIST_SNARL_PATHS

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <functional>
#include <iostream>
#include <chrono>
#include <cassert>
#include <regex>
#include <stdexcept>
#include <utility>
#include <bdsg/internal/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>

using namespace std;
using namespace bdsg;

using Snarl = bdsg::SnarlDistanceIndex::NetHandle;

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
    void addNodeHandle(const handle_t& node_h, const BaseHandleGraph& stree);

    // Get the string representation of the path
    std::string print() const;

    // Flip the path orientation
    void flip();

    // Get the size of the path
    size_t size() const;

    // Count the number of reversed nodes
    size_t nreversed() const;
};

// Function to split paths using regex
vector<string> split_paths(const string& path);

// Function to get length of a node
int length_node(const PathHandleGraph& pg, int node_id);

// Function to calculate the type of variant
vector<string> calcul_type_variant(const vector<vector<int>>& list_list_length_paths);

// Function to check threshold
void check_threshold(double proportion);

// Function to find snarl ID
string find_snarl_id(const SnarlTree& stree, const Snarl& snarl);

// Function to follow edges
void follow_edges(
    const SnarlTree& stree,
    vector<vector<string>>& finished_paths,
    const vector<string>& path,
    vector<vector<string>>& paths,
    const PathHandleGraph& pg
);

// Function to save snarls
vector<Snarl> save_snarls(SnarlDistanceIndex& stree, const Snarl& root);

// Function to parse the graph and tree from files
void parse_graph_tree(const string& pg_file, const string& dist_file, 
                      PackedGraph& pg, SnarlDistanceIndex& stree, Snarl& root);

// Function to fill pretty paths
pair<vector<string>, vector<vector<string>>> fill_pretty_paths(SnarlDistanceIndex& stree, 
                                                               PackedGraph& pg, 
                                                               const vector<vector<Snarl>>& finished_paths);

// Function to write the header to the output file
void write_header_output(const string& output_file);

// Function to write output to the file
void write_output(const string& output_file, const string& snarl_id, 
                  const vector<string>& pretty_paths, const vector<string>& type_variants);

// Function to loop over snarls and write output
void loop_over_snarls_write(SnarlDistanceIndex& stree, 
                            const vector<Snarl>& snarls, 
                            PackedGraph& pg, 
                            const string& output_file, 
                            const string& output_snarl_not_analyse, 
                            int time_threshold = 10);

#endif 
