#include "list_snarl_paths.hpp"

using namespace std;
using namespace bdsg;
using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

Path::Path() {}

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
            node_s.replace(pos, 25, "");
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
            size_t pos_2;
            while ((pos_2 = str.find(sub)) != std::string::npos) {
                str.erase(pos_2, sub.length());
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
pair<vector<string>, size_t> calcul_pos_type_variant(const vector<vector<string>>& list_list_length_paths) {
    vector<string> list_type_variant;
    size_t padding = 0;
    bool just_snp = true;

    for (const auto& path_lengths : list_list_length_paths) {
        if (path_lengths.size() > 3 || path_lengths[1] == "_") { // Case snarl in snarl / Indel
            list_type_variant.push_back("CPX"); // COMPLEX
            just_snp = false;
        } else if (path_lengths.size() == 3) { // Case simple path len 3
            if (path_lengths[1].size() == 1) {
                list_type_variant.push_back(path_lengths[1]); // add node str snp 
            } else {
                // vector string path_lengths
                string ins_seq = (path_lengths[1].size() > 3) ? "INS" : path_lengths[1];
                list_type_variant.push_back(ins_seq);
                just_snp = false;
            }
        } else if (path_lengths.size() == 2) { // Deletion
            list_type_variant.push_back("DEL");
            just_snp = false;
        } else { // Case path_lengths is empty
            cerr << "path_lengths is empty" << endl;
        }
    }

    // add +1 in pos for just SNP present in snarl
    if (just_snp) {padding = 1;}

    return {list_type_variant, padding};
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

std::tuple<std::unique_ptr<bdsg::SnarlDistanceIndex>, 
           std::unique_ptr<bdsg::PackedGraph>, 
           handlegraph::net_handle_t, 
           std::unique_ptr<bdsg::PackedPositionOverlay>>
parse_graph_tree(const std::string& pg_file, const std::string& dist_file) {
    
    // Load graph
    auto pg = std::make_unique<bdsg::PackedGraph>();
    pg->deserialize(pg_file);

    // Load snarl tree
    auto stree = std::make_unique<bdsg::SnarlDistanceIndex>();
    stree->deserialize(dist_file);

    // PackedPositionOverlay takes a pointer to pg
    auto pp_overlay = std::make_unique<bdsg::PackedPositionOverlay>(pg.get());

    // Get root of snarl tree
    handlegraph::net_handle_t root = stree->get_root();

    return std::make_tuple(std::move(stree), std::move(pg), root, std::move(pp_overlay));
}

void follow_edges(SnarlDistanceIndex& stree, 
                vector<vector<net_handle_t>>& finished_paths,
                vector<net_handle_t>& path,
                vector<vector<net_handle_t>>& paths, 
                PackedGraph& pg) {
  
    auto add_to_path = [&](const net_handle_t& next_child) {
        if (stree.is_sentinel(next_child)) {
            // If this is the bound of the snarl then we're done
            finished_paths.emplace_back(path);
            finished_paths.back().push_back(next_child);
        } else {
            for (const auto& i : path) {
                // Case where we find a loop
                if (stree.net_handle_as_string(i) == stree.net_handle_as_string(next_child)) {
                    return false;
                }
            }
            paths.emplace_back(path);
            paths.back().push_back(next_child);
        }
        return true;
    };

    // Follow edges from the last element in path
    if (!path.empty()) {
        stree.follow_net_edges(path.back(), &pg, false, add_to_path);
    }
}

vector<tuple<net_handle_t, string, size_t>> save_snarls(
                                SnarlDistanceIndex& stree, 
                                net_handle_t& root,
                                PackedGraph& pg, 
                                unordered_set<string>& ref_chr,
                                PackedPositionOverlay& ppo) {

    vector<tuple<net_handle_t, string, size_t>> snarls;
    unordered_map<string, pair<string, size_t>> snarls_pos;

    // Given a node handle (dist index), return a position on a chr reference path
    auto get_node_position = [&](net_handle_t node) -> pair<string, size_t> { // node : net_handle_t
        handle_t node_h = stree.get_handle(node, &pg);
        pair<string, size_t> ret_pos; // pair<string, size_t> path_name, position

        auto step_callback = [&](const step_handle_t& step_handle) {
            path_handle_t path_handle = pg.get_path_handle_of_step(step_handle);
            string chr_path = pg.get_path_name(path_handle);
            
            // check if chr_path is in ref_chr
            if (ref_chr.find(chr_path) != ref_chr.end()) {
                ret_pos.first = chr_path;
                ret_pos.second = ppo.get_position_of_step(step_handle) + stree.node_length(node); // position + length_node
                // cout << stree.net_handle_as_string(node) << " : " << pg.get_sequence(stree.get_handle(node, &pg)) << endl;
                return (false); // Stop iteration once a reference chr is found
            }
            return (true); // Continue iteration
        };

        pg.for_each_step_on_handle(node_h, step_callback);
        return ret_pos;
    };

    auto get_net_start_position = [&](net_handle_t net) -> pair<string, size_t> {

        if (stree.is_node(net)) {
            return get_node_position(net);
        }

        net_handle_t bnode1 = stree.get_bound(net, true, false);
        pair<string, size_t> bnode1_p = get_node_position(bnode1);
        net_handle_t bnode2 = stree.get_bound(net, false, false);
        pair<string, size_t> bnode2_p = get_node_position(bnode2);

        // Check if the string part of the pair is empty
        if (bnode1_p.first.empty()) return bnode1_p;
        if (bnode2_p.first.empty()) return bnode2_p;

        assert(bnode1_p.first == bnode2_p.first); // Ensure they are on the same reference path

        // Return the boundary with the smaller numerical position
        return (bnode1_p.second < bnode2_p.second) ? bnode1_p : bnode2_p;
    };

    function<void(net_handle_t)> save_snarl_tree_node;
    save_snarl_tree_node = [&](net_handle_t net) {
        pair<string, size_t> snarl_pos = get_net_start_position(net);
        if (snarl_pos.first.empty()) {
            auto par_net = stree.get_parent(net);
            snarl_pos = snarls_pos[stree.net_handle_as_string(par_net)];
        }
 
        snarls_pos[stree.net_handle_as_string(net)] = snarl_pos;
        if (stree.is_snarl(net)) {
            snarls.push_back(std::make_tuple(net, snarl_pos.first, snarl_pos.second));
        }

        if (!stree.is_node(net) && !stree.is_sentinel(net)) {
            stree.for_each_child(net, save_snarl_tree_node);
        }
    };
    
    stree.for_each_child(root, save_snarl_tree_node);
    cout << "Number of snarls : " << snarls.size() << endl;

    return snarls;
}

tuple<vector<string>, vector<string>, size_t> fill_pretty_paths(
    SnarlDistanceIndex& stree, 
    PackedGraph& pg, 
    vector<vector<net_handle_t>>& finished_paths) {
    
    vector<string> pretty_paths;
    vector<vector<string>> seq_net_paths;

    for (const auto& path : finished_paths) {
        Path ppath;
        vector<string> seq_net;

        for (auto net : path) {
            if (stree.is_sentinel(net)) {
                net = stree.get_node_from_sentinel(net);
            }
            
            // Node case
            if (stree.is_node(net)) {
                ppath.addNodeHandle(net, stree);
                //length_net.push_back(stree.node_length(net));
                
                nid_t node_start_id = stree.node_id(net);
                handle_t node_handle = pg.get_handle(node_start_id);
                string seq_node = pg.get_sequence(node_handle);
                seq_net.push_back(seq_node);
            }

            // Trivial chain case
            else if (stree.is_trivial_chain(net)) {
                ppath.addNodeHandle(net, stree);
                auto stn_start = stree.starts_at_start(net) ? stree.get_bound(net, false, true) : stree.get_bound(net, true, true);
                auto node_start_id = stree.node_id(stn_start);
                auto net_trivial_chain = pg.get_handle(node_start_id);
                string seq_trivial_chain = pg.get_sequence(net_trivial_chain);
                seq_net.push_back(seq_trivial_chain);
            }

            // Chain case
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
                seq_net.push_back("_");
            }
        }

        if (ppath.nreversed() > ppath.size() / 2) {
            ppath.flip();
            std::reverse(seq_net.begin(), seq_net.end());
        }
        pretty_paths.push_back(ppath.print());
        seq_net_paths.push_back(seq_net);
    }

    // pair<vector<string>, size_t>
    auto [type_variants, length_first_variant] = calcul_pos_type_variant(seq_net_paths);
    return std::make_tuple(pretty_paths, type_variants, length_first_variant);
}

// {chr : matrix(snarl, paths, chr, pos, type)}
std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> loop_over_snarls_write(
        SnarlDistanceIndex& stree,
        vector<tuple<net_handle_t, string, size_t>>& snarls,
        PackedGraph& pg, 
        const string& output_file,
        const string& output_snarl_not_analyse,
        size_t children_threshold = 50, 
        bool bool_return = true) {

    ofstream out_snarl(output_file);
    ofstream out_fail(output_snarl_not_analyse);
    
    out_snarl << "chr\tpos\tsnarl\tpaths\ttype\n";
    out_fail << "snarl\treason\n";
    
    std::vector<std::tuple<string, vector<string>, string, vector<string>>> snarl_paths;
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> chr_snarl_matrix;
    size_t paths_number_analysis = 0;
    size_t itr_thresold = 10000; // TODO change it 
    string save_chr = "";

    std::vector<size_t> children = {0};
    auto count_children = [&](net_handle_t net) {
        children[0] += 1;
        return true;
    };

    for (const auto& snarl_path_pos : snarls) {
        net_handle_t snarl = std::get<0>(snarl_path_pos);
        size_t itr = 0;
        string snarl_id = find_snarl_id(stree, snarl);
        bool not_break = true;
        children = {0}; // re-initialise the children vec
        
        stree.for_each_child(snarl, count_children);
        if (children[0] > children_threshold) {
            out_fail << snarl_id << "\ttoo_many_children\n";
            continue;
        }
        
        vector<vector<net_handle_t>> paths = {{stree.get_bound(snarl, false, true)}};
        vector<vector<net_handle_t>> finished_paths;

        while (!paths.empty()) {
            auto path = paths.back();
            paths.pop_back();

            if (itr > itr_thresold) {
                out_fail << snarl_id << "\titeration_calculation_out\n";
                not_break = false;
                break;
            }
            follow_edges(stree, finished_paths, path, paths, pg);
            itr++;
        }

        if (not_break) {
            // pair<vector<string>, vector<string>>
            auto [pretty_paths, type_variants, padding] = fill_pretty_paths(stree, pg, finished_paths);
            std::ostringstream pretty_paths_stream, type_variants_stream;

            // Convert pretty_paths (vector<string>) into a comma-separated string
            for (size_t i = 0; i < pretty_paths.size(); ++i) {
                if (i > 0) pretty_paths_stream << ",";
                pretty_paths_stream << pretty_paths[i];
            }
 
            for (size_t i = 0; i < type_variants.size(); ++i) {
                if (i > 0) type_variants_stream << ",";  // Add a comma and space between strings
                type_variants_stream << type_variants[i];
            }

            // chromosome   position    snarl_id    paths    type
            string chr = std::get<1>(snarl_path_pos);
            string pos = std::to_string(std::get<2>(snarl_path_pos)+padding);
            paths_number_analysis += pretty_paths.size();

            if (bool_return) {
                out_snarl << chr << "\t" << pos
                    << "\t" << snarl_id << "\t" << pretty_paths_stream.str() 
                    << "\t" << type_variants_stream.str() << "\n";
            } else {
                // case new chr
                if (chr != save_chr && !save_chr.empty()) {
                    cout << "cleaning" << endl; 
                    chr_snarl_matrix[save_chr] = std::move(snarl_paths);
                    snarl_paths.clear();
                }
                save_chr = chr;
                snarl_paths.push_back(std::make_tuple(snarl_id, pretty_paths, pos, type_variants));
            }
        }
    }

    // last chr adding, but only if save_chr is not empty
    if (!save_chr.empty()) {
        chr_snarl_matrix[save_chr] = std::move(snarl_paths);
    }

    // Print the size of snarl_paths
    cout << "Final snarl_paths size : " << paths_number_analysis << endl;

    // Print chr_snarl_matrix
    for (const auto& chr_snarl : chr_snarl_matrix) {
        cout << "chr : " << chr_snarl.first << ", number of snarl : " << chr_snarl.second.size() << endl;
    }

    return {chr_snarl_matrix};
}
