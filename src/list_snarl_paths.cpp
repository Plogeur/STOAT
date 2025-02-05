#include "list_snarl_paths.hpp"

using namespace std;
using namespace bdsg;

using Snarl = bdsg::SnarlDistanceIndex::NetHandle;

// Add a node with known orientation
void addNode(const std::string& node, char orient) {
    nodes.push_back(node);
    orients.push_back(orient);
}

// Add a node handle and extract information using the string representation
void addNodeHandle(const handle_t& node_h, const BaseHandleGraph& stree) {
    std::string node_s = stree.net_handle_as_string(node_h);

    // Handle trivial chain modifications
    if (stree.is_trivial_chain(node_h)) {
        size_t pos;
        while ((pos = node_s.find(" pretending to be a chain")) != std::string::npos) {
            node_s.replace(pos, 23, "");
        }
        while ((pos = node_s.find(" in a simple snarl")) != std::string::npos) {
            node_s.replace(pos, 19, "");
        }
    }

    // Parse node info
    size_t pos = node_s.find("node ");
    if (pos != std::string::npos) {
        node_s.erase(pos, 5);
    }

    char node_o = '>';
    if (node_s.find("rev") != std::string::npos) {
        node_o = '<';
    }

    auto removeSubstrings = [](std::string& str, const std::vector<std::string>& substrings) {
        for (const auto& sub : substrings) {
            size_t pos;
            while ((pos = str.find(sub)) != std::string::npos) {
                str.erase(pos, sub.length());
            }
        }
    };

    removeSubstrings(node_s, {"rev", "fd"});

    // Add node to path
    nodes.push_back(node_s);
    orients.push_back(node_o);
}

// Get the string representation of the path
std::string print() const {
    std::string out_path;
    for (size_t i = 0; i < nodes.size(); ++i) {
        out_path += orients[i] + nodes[i];
    }
    return out_path;
}

// Flip the path orientation
void flip() {
    std::reverse(nodes.begin(), nodes.end());
    std::reverse(orients.begin(), orients.end());
    for (size_t i = 0; i < orients.size(); ++i) {
        if (nodes[i] == "*") {
            continue;
        }
        orients[i] = (orients[i] == '>') ? '<' : '>';
    }
}

// Get the size of the path
size_t size() const {
    return nodes.size();
}

// Count the number of reversed nodes
size_t nreversed() const {
    return std::count(orients.begin(), orients.end(), '<');
}

// Function to split paths using regex
vector<string> split_paths(const string& path) {
    regex re("\\d+");
    sregex_iterator begin(path.begin(), path.end(), re), end;
    vector<string> result;
    for (auto it = begin; it != end; ++it) {
        result.push_back(it->str());
    }
    return result;
}

// Function to get length of a node
int length_node(const PathHandleGraph& pg, int node_id) {
    return pg.get_length(node_id); // Assuming `get_length` exists in bdsg
}

// Function to calculate the type of variant
vector<string> calcul_type_variant(const vector<vector<int>>& list_list_length_paths) {
    vector<string> list_type_variant;

    for (const auto& path_lengths : list_list_length_paths) {
        if (path_lengths.size() > 3 || path_lengths[1] == -1) { // Case snarl in snarl / Indel
            list_type_variant.push_back("COMPLEX");
        } else if (path_lengths.size() == 3) { // Case simple path len 3
            list_type_variant.push_back((path_lengths[1] == 1) ? "SNP" : "INS");
        } else { // Deletion
            list_type_variant.push_back("DEL");
        }
    }

    return list_type_variant;
}

// Function to check threshold
void check_threshold(double proportion) {
    if (proportion <= 0) {
        throw invalid_argument("Proportion value must be >0.");
    }
}

// Function to find snarl ID
string find_snarl_id(const SnarlTree& stree, const Snarl& snarl) {
    auto sstart = stree.get_bound(snarl, false, true);
    sstart = stree.get_node_from_sentinel(sstart);
    auto send = stree.get_bound(snarl, true, true);
    send = stree.get_node_from_sentinel(send);

    return to_string(stree.node_id(send)) + "_" + to_string(stree.node_id(sstart));
}

// Function to follow edges
void follow_edges(
    const SnarlTree& stree,
    vector<vector<string>>& finished_paths,
    const vector<string>& path,
    vector<vector<string>>& paths,
    const PathHandleGraph& pg
) {
    auto add_to_path = [&](const auto& next_child) -> bool {
        if (stree.is_sentinel(next_child)) {
            // If this is the bound of the snarl, we're done
            finished_paths.emplace_back(path);
            finished_paths.back().push_back(stree.net_handle_as_string(next_child));
        } else {
            for (const auto& i : path) {
                // Case where we find a loop
                if (stree.net_handle_as_string(i) == stree.net_handle_as_string(next_child)) {
                    return false;
                }
            }
            paths.emplace_back(path);
            paths.back().push_back(stree.net_handle_as_string(next_child));
        }
        return true;
    };

    // Follow the net edges from the last element in the path
    stree.follow_net_edges(stree.net_handle_from_string(path.back()), pg, false, add_to_path);
}

vector<vector<string>> save_snarls(SnarlTree& stree, SnarlTree::NetHandle root,
                                   PackedGraph& pg, unordered_set<string>& ref_paths,
                                   PathPositionOverlay& ppo) {
    vector<vector<string>> snarls;
    unordered_map<string, Position> snarls_pos;

    // Given a node handle (dist index) return a position on a reference path
    auto get_node_position = [&](typename SnarlTree::NetHandle node) -> Position {
        auto node_h = stree.get_handle(node, pg);
        Position ret_pos;

        pg.for_each_step_on_handle(node_h, [&](step_handle_t step_handle) {
            path_handle_t path_handle = pg.get_path_handle_of_step(step_handle);
            string path_name = pg.get_path_name(path_handle);
            if (ref_paths.count(path_name)) {
                size_t position = ppo.get_position_of_step(step_handle);
                ret_pos.push_back(path_name);
                ret_pos.push_back(to_string(position));
                return false; // Stop iteration
            }
            return true;
        });
        return ret_pos;
    };

    auto get_net_start_position = [&](typename SnarlTree::NetHandle net) -> Position {
        if (stree.is_node(net)) {
            return get_node_position(net);
        }
        auto bnode1 = stree.get_bound(net, true, false);
        auto bnode1_p = get_node_position(bnode1);
        auto bnode2 = stree.get_bound(net, false, false);
        auto bnode2_p = get_node_position(bnode2);

        if (bnode1_p.empty()) return bnode1_p;
        if (bnode2_p.empty()) return bnode2_p;
        assert(bnode1_p[0] == bnode2_p[0]); // Ensure they are on the same ref path
        return (stoi(bnode1_p[1]) < stoi(bnode2_p[1])) ? bnode1_p : bnode2_p;
    };

    function<void(typename SnarlTree::NetHandle)> save_snarl_tree_node;
    save_snarl_tree_node = [&](typename SnarlTree::NetHandle net) {
        Position snarl_pos = get_net_start_position(net);
        if (snarl_pos.empty()) {
            auto par_net = stree.get_parent(net);
            snarl_pos = snarls_pos[stree.net_handle_as_string(par_net)];
        }
        snarls_pos[stree.net_handle_as_string(net)] = snarl_pos;
        if (stree.is_snarl(net)) {
            snarls.push_back({stree.net_handle_as_string(net), snarl_pos[0], snarl_pos[1]});
        }
        if (!stree.is_node(net) && !stree.is_sentinel(net)) {
            stree.for_each_child(net, save_snarl_tree_node);
        }
    };

    stree.for_each_child(root, save_snarl_tree_node);
    return snarls;
}
struct GraphTree {
    SnarlDistanceIndex stree;
    PackedGraph pg;
    SnarlDistanceIndex::NetHandle root;
    PackedPositionOverlay pp_overlay;
};

GraphTree parse_graph_tree(const string& pg_file, const string& dist_file) {
    // Load graph and snarl tree
    GraphTree graph_tree;
    
    graph_tree.pg.deserialize(pg_file);
    graph_tree.stree.deserialize(dist_file);
    graph_tree.pp_overlay = PackedPositionOverlay(graph_tree.pg);

    // Get root snarl
    graph_tree.root = graph_tree.stree.get_root();
    
    return graph_tree;
}

pair<vector<string>, vector<vector<string>>> fill_pretty_paths(
    SnarlTree& stree, PackedGraph& pg, 
    vector<vector<typename SnarlTree::NetHandle>> finished_paths) {
    
    vector<string> pretty_paths;
    vector<vector<string>> length_net_paths;

    for (const auto& path : finished_paths) {
        Path ppath;
        vector<string> length_net;

        for (auto net : path) {
            if (stree.is_sentinel(net)) {
                net = stree.get_node_from_sentinel(net);
            }

            if (stree.is_node(net)) {
                ppath.addNodeHandle(net, stree);
                length_net.push_back(to_string(stree.node_length(net)));
            }
            else if (stree.is_trivial_chain(net)) {
                ppath.addNodeHandle(net, stree);
                auto stn_start = stree.get_bound(net, false, true);
                auto node_start_id = stree.node_id(stn_start);
                auto net_trivial_chain = pg.get_handle(node_start_id);
                length_net.push_back(to_string(pg.get_length(net_trivial_chain)));
            }
            else if (stree.is_chain(net)) {
                typename SnarlTree::NetHandle nodl, nodr;
                if (stree.starts_at_start(net)) {
                    nodl = stree.get_bound(net, false, true);
                    nodr = stree.get_bound(net, true, false);
                } else {
                    nodl = stree.get_bound(net, true, true);
                    nodr = stree.get_bound(net, false, false);
                }
                ppath.addNodeHandle(nodl, stree);
                ppath.addNode('*', '>');
                ppath.addNodeHandle(nodr, stree);
                length_net.push_back("-1");
            }
        }

        if (ppath.nreversed() > ppath.size() / 2) {
            ppath.flip();
        }
        pretty_paths.push_back(ppath.print());
        length_net_paths.push_back(length_net);
    }

    return {pretty_paths, length_net_paths};
}

void write_header_output(const string& output_file) {
    ofstream outf(output_file);
    outf << "snarl\tpaths\ttype\n";
}

void write_output(const string& output_file, const string& snarl_id, 
                  const vector<string>& pretty_paths, const vector<string>& type_variants) {
    ofstream outf(output_file, ios::app);
    outf << snarl_id << "\t";
    for (size_t i = 0; i < pretty_paths.size(); ++i) {
        outf << pretty_paths[i] << (i + 1 < pretty_paths.size() ? "," : "");
    }
    outf << "\t";
    for (size_t i = 0; i < type_variants.size(); ++i) {
        outf << type_variants[i] << (i + 1 < type_variants.size() ? "," : "");
    }
    outf << "\n";
}

void loop_over_snarls_write(SnarlTree& stree, vector<vector<typename SnarlTree::NetHandle>> snarls,
                        PackedGraph& pg, const string& output_file,
                        const string& output_snarl_not_analyse,
                        int children_threshold = 50, bool bool_return = true) {
    ofstream out_snarl(output_file);
    ofstream out_fail(output_snarl_not_analyse);
    
    out_snarl << "snarl\tpaths\ttype\tchr\tpos\n";
    out_fail << "snarl\treason\n";
    chrono::seconds time_threshold(2);
    
    for (const auto& snarl_path_pos : snarls) {
        auto snarl = snarl_path_pos[0];
        auto snarl_time = chrono::steady_clock::now();
        string snarl_id = find_snarl_id(stree, snarl);
        bool not_break = true;
        int children = 0;
        
        auto count_children = [&](auto net) {
            children++;
            return true;
        };
        
        stree.for_each_child(snarl, count_children);
        if (children > children_threshold) {
            out_fail << snarl_id << "\ttoo_many_children\n";
            continue;
        }
        
        vector<vector<typename SnarlTree::NetHandle>> paths = {{stree.get_bound(snarl, false, true)}};
        vector<vector<typename SnarlTree::NetHandle>> finished_paths;
        
        while (!paths.empty()) {
            auto path = paths.back();
            paths.pop_back();
            
            if (chrono::steady_clock::now() - snarl_time > time_threshold) {
                out_fail << snarl_id << "\ttime_calculation_out\n";
                not_break = false;
                break;
            }
            
            follow_edges(stree, finished_paths, path, paths, pg);
        }
        
        if (not_break) {
            auto [pretty_paths, type_variants] = fill_pretty_paths(stree, pg, finished_paths);
            
            out_snarl << snarl_id << "\t" << join(pretty_paths, ",") << "\t"
                      << join(type_variants, ",") << "\t" << snarl_path_pos[1]
                      << "\t" << snarl_path_pos[2] << "\n";
        }
    }
}