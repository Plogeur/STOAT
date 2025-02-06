#include "list_snarl_paths.hpp"

using namespace std;
using namespace bdsg;
using bdsg::HandleGraph;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

// Add a node with known orientation
void Path::addNode(const std::string& node, char orient) {
    nodes.push_back(node);
    orients.push_back(orient);
}

// Add a node handle and extract information using the string representation
void Path::addNodeHandle(const net_handle_t& node_h, const SnarlDistanceIndex& stree) {
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
std::string Path::print() const {
    std::string out_path;
    for (size_t i = 0; i < nodes.size(); ++i) {
        out_path += orients[i] + nodes[i];
    }
    return out_path;
}

// Flip the path orientation
void Path::flip() {
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
size_t Path::size() const {
    return nodes.size();
}

// Count the number of reversed nodes
size_t Path::nreversed() const {
    return std::count(orients.begin(), orients.end(), '<');
}

// Function to calculate the type of variant
vector<string> calcul_type_variant(const vector<vector<string>>& list_list_length_paths) {
    vector<string> list_type_variant;

    for (const auto& path_lengths : list_list_length_paths) {
        if (path_lengths.size() > 3 || path_lengths[1] == "-1") { // Case snarl in snarl / Indel
            list_type_variant.push_back("COMPLEX");
        } else if (path_lengths.size() == 3) { // Case simple path len 3
            list_type_variant.push_back((path_lengths[1] == "1") ? "SNP" : "INS");
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

string find_snarl_id(SnarlDistanceIndex& stree, net_handle_t& snarl) {
    // Get start and end boundary nodes for the snarl
    auto sstart = stree.get_bound(snarl, false, true);  // False for the left boundary
    auto send = stree.get_bound(snarl, true, true);     // True for the right boundary

    // Convert the sentinels into nodes
    auto start_node = stree.get_node_from_sentinel(sstart);
    auto end_node = stree.get_node_from_sentinel(send);

    // Get the node IDs from SnarlDistanceIndex
    auto start_node_id = stree.node_id(start_node);
    auto end_node_id = stree.node_id(end_node);

    // Construct the snarl ID as "end_node_id_start_node_id"
    std::stringstream snarl_id;
    snarl_id << end_node_id << "_" << start_node_id;

    return snarl_id.str();  // Return the generated snarl ID as a string
}

void follow_edges(SnarlDistanceIndex& stree, vector<vector<net_handle_t>>& finished_paths,
                  vector<net_handle_t>& path,
                  vector<vector<net_handle_t>>& paths, PackedGraph& pg) {
}

pair<vector<tuple<string, vector<string>, string, string, string>>, int> loop_over_snarls_write(
        SnarlDistanceIndex& stree, 
        vector<tuple<net_handle_t, string, size_t>> snarls,
        PackedGraph& pg, const string& output_file,
        const string& output_snarl_not_analyse,
        int children_threshold = 50, bool bool_return = true) {

    ofstream out_snarl(output_file);
    ofstream out_fail(output_snarl_not_analyse);
    
    out_snarl << "snarl\tpaths\ttype\tchr\tpos\n";
    out_fail << "snarl\treason\n";
    
    vector<tuple<string, vector<string>, string, string, string>> snarl_paths;
    int paths_number_analysis = 0;
    chrono::seconds time_threshold(2);
    
    for (const auto& snarl_path_pos : snarls) {
        net_handle_t snarl = std::get<0>(snarl_path_pos);
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
        
        vector<vector<net_handle_t>> paths = {{stree.get_bound(snarl, false, true)}};
        vector<vector<net_handle_t>> finished_paths;

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
            // pair<vector<string>, vector<string>>
            auto [pretty_paths, type_variants] = fill_pretty_paths(stree, pg, finished_paths);
            std::ostringstream pretty_paths_stream, type_variants_stream;

            // Convert pretty_paths (vector<string>) into a comma-separated string
            for (size_t i = 0; i < pretty_paths.size(); ++i) {
                if (i > 0) pretty_paths_stream << ",";
                pretty_paths_stream << pretty_paths[i];
            }

            // Convert type_variants to a comma-separated string
            for (size_t i = 0; i < type_variants.size(); ++i) {
                if (i > 0) type_variants_stream << ","; // Separate vectors with commas

                // Convert inner vector (type_variants[i]) to a comma-separated string
                for (size_t j = 0; j < type_variants[i].size(); ++j) {
                    if (j > 0) type_variants_stream << " "; // Space-separated words within a vector
                    type_variants_stream << type_variants[i][j];
                }
            }

            out_snarl << snarl_id << "\t" << pretty_paths_stream.str() << "\t"
                    << type_variants_stream.str() << "\t" << std::get<1>(snarl_path_pos)
                    << "\t" << std::to_string(std::get<2>(snarl_path_pos)) << "\n";

            if (bool_return) {
                snarl_paths.emplace_back(snarl_id, pretty_paths, type_variants_stream.str(),
                                        std::get<1>(snarl_path_pos), std::to_string(std::get<2>(snarl_path_pos)));
            }

            paths_number_analysis += pretty_paths.size();
        }
    }
    
    return {snarl_paths, paths_number_analysis};
}

vector<tuple<net_handle_t, string, size_t>> save_snarls(SnarlDistanceIndex& stree, const net_handle_t& root,
                                 PackedGraph& pg, unordered_set<string>& ref_paths,
                                 PackedPositionOverlay& ppo) {
    vector<tuple<net_handle_t, string, size_t>> snarls;
    unordered_map<string, vector<string>> snarls_pos;

    // Given a node handle (dist index), return a position on a reference path
    auto get_node_position = [&](auto node) -> vector<string> { // node : net_handle_t
        handle_t node_h = stree.get_handle(node, pg);
        vector<string> ret_pos(2); // pair<string, size_t>

        pg.for_each_step_on_handle(node_h, [&](auto step_handle) {
            path_handle_t path_handle = pg.get_path_handle_of_step(step_handle);
            string path_name = pg.get_path_name(path_handle);

            if (ref_paths.count(path_name)) {
                size_t position = ppo.get_position_of_step(step_handle);
                ret_pos[0] = path_name;
                ret_pos[1] = position;
                return false; // Stop iteration once a reference path is found
            }
            return true; // Continue iteration
        });
        return ret_pos;
    };

    auto get_net_start_position = [&](net_handle_t net) -> vector<string> {
        if (stree.is_node(net)) {
            return get_node_position(net);
        }

        net_handle_t bnode1 = stree.get_bound(net, true, false);
        vector<string> bnode1_p = get_node_position(bnode1);
        net_handle_t bnode2 = stree.get_bound(net, false, false);
        vector<string> bnode2_p = get_node_position(bnode2);

        // If one of the boundaries is not on a reference path, return it
        if (bnode1_p.empty()) return bnode1_p;
        if (bnode2_p.empty()) return bnode2_p;

        assert(bnode1_p[0] == bnode2_p[0]); // Ensure they are on the same reference path

        // Return the boundary with the smaller numerical position
        return (std::stoi(bnode1_p[1]) < std::stoi(bnode2_p[1])) ? bnode1_p : bnode2_p;
    };

    function<void(net_handle_t)> save_snarl_tree_node;
    save_snarl_tree_node = [&](net_handle_t net) {
        vector<string> snarl_pos = get_net_start_position(net);
        if (snarl_pos.empty()) {
            auto par_net = stree.get_parent(net);
            snarl_pos = snarls_pos[stree.net_handle_as_string(par_net)];
        }

        snarls_pos[stree.net_handle_as_string(net)] = snarl_pos;
        if (stree.is_snarl(net)) {
            snarls.push_back(std::make_tuple(net, snarl_pos[0], std::stoull(snarl_pos[1])));
        }

        if (!stree.is_node(net) && !stree.is_sentinel(net)) {
            stree.for_each_child(net, save_snarl_tree_node);
        }
    };
    
    stree.for_each_child(root, save_snarl_tree_node);
    return snarls;
}

std::tuple<bdsg::PackedGraph, bdsg::SnarlDistanceIndex, bdsg::PackedPositionOverlay, net_handle_t>
parse_graph_tree(const std::string& pg_file, const std::string& dist_file) {
    // Load graph and snarl tree
    bdsg::PackedGraph pg;
    pg.deserialize(pg_file);
    
    bdsg::SnarlDistanceIndex stree;
    stree.deserialize(dist_file);

    bdsg::PathHandleGraph* path_graph = &pg;
    bdsg::PackedPositionOverlay pp_overlay(path_graph);
    // /Users/aliasmatis/Desktop/stoat_cxx/tests/src/list_snarl_paths.cpp:302:33: error: no matching constructor for initialization of 'bdsg::PackedPositionOverlay'
    // bdsg::PackedPositionOverlay pp_overlay(pg);

    // Get the root of the snarl tree
    auto root = stree.get_root();
    
    return std::make_tuple(pg, stree, pp_overlay, root);
}

pair<vector<string>, vector<string>> fill_pretty_paths(
    SnarlDistanceIndex& stree, PackedGraph& pg, 
    vector<vector<net_handle_t>>& finished_paths) {
    
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
                net_handle_t nodl, nodr;
                if (stree.starts_at_start(net)) {
                    nodl = stree.get_bound(net, false, true);
                    nodr = stree.get_bound(net, true, false);
                } else {
                    nodl = stree.get_bound(net, true, true);
                    nodr = stree.get_bound(net, false, false);
                }
                ppath.addNodeHandle(nodl, stree);
                ppath.addNode("*", '>');
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

    vector<string> type_variants = calcul_type_variant(length_net_paths);
    return {pretty_paths, type_variants};
}
