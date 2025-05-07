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
bool Path::addNodeHandle(const net_handle_t& node_h, const SnarlDistanceIndex& stree) {
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
    return node_o == '>' ? true : false;
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
// tuple<string, size_t, size_t, size_t>
// seq_net, minimum_distance, maximun_distance, size_path, sum_path
vector<string> calcul_pos_type_variant(const vector<tuple<string, size_t, size_t, size_t, size_t>>& list_length_paths) {
    vector<string> list_type_variant;

    for (const auto& tuple_info : list_length_paths) {
        size_t path_length = std::get<3>(tuple_info);
        size_t sum_path = std::get<4>(tuple_info);
        if (path_length > 3) {
            if (sum_path == 0) { // Case complex 
                string complex = to_string(std::get<1>(tuple_info)) + "/" + to_string(std::get<2>(tuple_info));
                list_type_variant.push_back(complex);
            } else { // Case multiple nodes (ex : INS+SNP+...)
                list_type_variant.push_back(std::to_string(sum_path));
            }

        } else if (path_length == 3) { // Case simple path len 3
            string seq = std::get<0>(tuple_info);
            size_t seq_length = seq.length();
            list_type_variant.push_back(to_string(seq_length));

        } else if (path_length == 2) { // case Deletion
            list_type_variant.push_back("0");
        } else { // Case path_lengths is empty
            cerr << "path_lengths is empty" << endl;
        }
    }

    return list_type_variant;
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
            parse_graph_tree(const std::string& pg_file, 
                const std::string& dist_file) {
                
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

        cout << "stree.net_handle_as_string(next_child) : " << stree.net_handle_as_string(next_child) << endl;
        if (stree.is_sentinel(next_child)) {
            // If this is the bound of the snarl then we're done
            finished_paths.emplace_back(path);
            finished_paths.back().push_back(next_child);
        } else {
            // Case where we find a loop
            for (const auto& i : path) {
                cout << "stree.net_handle_as_string(i) : " << stree.net_handle_as_string(i) << endl;
                if (stree.net_handle_as_string(i) == stree.net_handle_as_string(next_child)) {
                    cout << "loop found" << endl;
                    return false;
                }
            }
            paths.emplace_back(path);
            paths.back().push_back(next_child);
        }
        cout << endl;
        return true;
    };

    // Follow edges from the last element in path
    if (!path.empty()) {
        stree.follow_net_edges(path.back(), &pg, false, add_to_path);
    }
}

vector<tuple<net_handle_t, string, size_t, size_t, bool>> save_snarls(
                                SnarlDistanceIndex& stree, 
                                net_handle_t& root,
                                PackedGraph& pg, 
                                unordered_set<string>& ref_chr,
                                PackedPositionOverlay& ppo) {

    vector<tuple<net_handle_t, string, size_t, size_t, bool>> snarls;
    unordered_map<string, tuple<string, size_t, size_t>> snarls_pos;
    size_t save_end_pos_ref = 0;

    // Given a node handle (dist index), return a position if on chr reference path
    auto get_node_position = [&](net_handle_t node) -> tuple<string, size_t, size_t> { // node : net_handle_t
        handle_t node_h = stree.get_handle(node, &pg);

        // path_name, position
        tuple<string, size_t, size_t> ret_pos;

        auto step_callback = [&](const step_handle_t& step_handle) {
            path_handle_t path_handle = pg.get_path_handle_of_step(step_handle);
            string chr_path = pg.get_path_name(path_handle);

            // check if chr_path is in ref_chr
            if (ref_chr.find(chr_path) != ref_chr.end()) {
                std::get<0>(ret_pos) = chr_path;
                size_t pos = ppo.get_position_of_step(step_handle);
                std::get<1>(ret_pos) = pos + stree.node_length(node); // position + length_node
                std::get<2>(ret_pos) = pos+1 ; // default end position

                return (false); // Stop iteration once a reference chr is found
            }
            return (true); // Continue iteration
        };

        pg.for_each_step_on_handle(node_h, step_callback);
        return ret_pos;
    };

    auto get_net_start_position = [&](net_handle_t net) -> tuple<string, size_t, size_t> {

        if (stree.is_node(net)) {
            return get_node_position(net);
        }

        net_handle_t bnode1 = stree.get_bound(net, true, false);
        tuple<string, size_t, size_t> bnode1_p = get_node_position(bnode1);

        net_handle_t bnode2 = stree.get_bound(net, false, false); // verify false true ?
        tuple<string, size_t, size_t> bnode2_p = get_node_position(bnode2);

        // Check if the string part of the pair is empty
        if (std::get<0>(bnode1_p).empty()) return bnode1_p;
        if (std::get<0>(bnode2_p).empty()) return bnode2_p;

        assert(std::get<0>(bnode1_p) == std::get<0>(bnode2_p)); // Ensure they are on the same reference path

        size_t start;
        size_t end;
        // smaller boundary is the start
        // larger boundary is the end
        if (std::get<1>(bnode1_p) < std::get<1>(bnode2_p)) {
            start = std::get<1>(bnode1_p);
            end = std::get<2>(bnode2_p);
        } else {
            start = std::get<1>(bnode2_p);
            end = std::get<2>(bnode1_p);
        }

        // tuple<string, size_t, size_t> snarl_start_end;
        return make_tuple(std::get<0>(bnode1_p), start, end);
    };

    function<void(net_handle_t)> save_snarl_tree_node;
    save_snarl_tree_node = [&](net_handle_t net) {

        tuple<string, size_t, size_t> snarl_pos = get_net_start_position(net);
        bool bool_ref = true;

        // if we couldn't find a position, use the parent's that we should have
        // found and saved earlier
        if (std::get<0>(snarl_pos).empty()) {
            auto par_net = stree.get_parent(net);
            snarl_pos = snarls_pos[stree.net_handle_as_string(par_net)];
            bool_ref = false;
        }

        // save this position
        snarls_pos[stree.net_handle_as_string(net)] = snarl_pos;

        // save snarl
        if (stree.is_snarl(net)) {
            // net_handle_t snarl, chr_ref, pos, is_on_ref_bool
            snarls.push_back(std::make_tuple(net, std::get<0>(snarl_pos), std::get<1>(snarl_pos), std::get<2>(snarl_pos), bool_ref));
        }

        // explore children
        if (!stree.is_node(net) && !stree.is_sentinel(net)) {
            stree.for_each_child(net, save_snarl_tree_node);
        }
    };

    stree.for_each_child(root, save_snarl_tree_node);
    cout << "Number of snarls : " << snarls.size() << endl;
    return snarls;
}

tuple<vector<string>, vector<string>> fill_pretty_paths(
    SnarlDistanceIndex& stree, 
    PackedGraph& pg, 
    vector<vector<net_handle_t>>& finished_paths) {
    
    // list of paths
    vector<string> pretty_paths;

    // seq_net, minimum_distance, maximum_distance
    vector<tuple<string, size_t, size_t, size_t, size_t>> seq_net_paths;

    for (const auto& path : finished_paths) {
        Path ppath;
        string seq_net;
        bool is_complex = false;
        size_t sum_path = 0;
        size_t minimum_distance=0;
        size_t maximun_distance=0;
        std::vector<size_t> size_node;
        size_node.resize(path.size());

        for (int i=0; i<path.size(); i++) {
            net_handle_t net = path[i];

            if (stree.is_sentinel(net)) {
                net = stree.get_node_from_sentinel(net);
            }

            // Node case
            if (stree.is_node(net)) {
                bool rev = ppath.addNodeHandle(net, stree);
                nid_t node_start_id = stree.node_id(net);
                handle_t node_handle = pg.get_handle(node_start_id);
                size_node[i] = pg.get_length(node_handle);
                if (ppath.size() == 2) { // add only the node seq in position 2 on the snarl (ex : X>P>Q, P is in position 2)
                    seq_net = pg.get_sequence(node_handle);
                }
            }

            // Trivial chain case
            else if (stree.is_trivial_chain(net)) {
                bool rev = ppath.addNodeHandle(net, stree);
                auto stn_start = stree.starts_at_start(net) ? stree.get_bound(net, false, true) : stree.get_bound(net, true, true);
                nid_t node_start_id = stree.node_id(stn_start);
                handle_t net_trivial_chain = pg.get_handle(node_start_id);
                size_node[i] = pg.get_length(net_trivial_chain);
                if (ppath.size() == 2) {
                    seq_net = pg.get_sequence(net_trivial_chain);
                }
            }

            // Chain case aka complexe
            else if (stree.is_chain(net)) {
                net_handle_t nodl, nodr;
                bool boundl;
                if (stree.starts_at_start(net)) {
                    cout << "start at start 1" << endl;
                    boundl = false;
                    nodl = stree.get_bound(net, false, true);
                    nodr = stree.get_bound(net, true, false);
                } else {
                    cout << "start at start 2" << endl;
                    boundl = true;
                    nodl = stree.get_bound(net, true, true);
                    nodr = stree.get_bound(net, false, false);
                }

                ppath.addNodeHandle(nodl, stree);
                ppath.addNode("*", '>');
                ppath.addNodeHandle(nodr, stree);

                // Get the size of the chain and return the distance (minimum and maximum)
                size_t complex_start_id = stree.node_id(nodl);
                size_t size_start_node = pg.get_length(pg.get_handle(complex_start_id));
                size_t complex_end_id = stree.node_id(nodr);
                size_t size_end_node = pg.get_length(pg.get_handle(complex_end_id));
                size_t size_chain = size_start_node + size_end_node;

                // boundl = true or false 
                size_t min_dist = stree.minimum_distance(complex_start_id, boundl, size_start_node, complex_end_id, boundl, 0);
                size_t max_dist = stree.maximum_distance(complex_start_id, boundl, size_start_node, complex_end_id, boundl, 0);
                minimum_distance = size_chain + min_dist;
                maximun_distance = size_chain + max_dist;
                // cout << "stree.node_id(nodl) : " << stree.node_id(nodl) << endl;
                // cout << "stree.node_id(nodr) : " << stree.node_id(nodr) << endl;
                // cout << "size_end_node = " << size_end_node << endl;
                // cout << "size_start_node = " << size_start_node << endl;
                // cout << "size_chain = " << size_chain << endl;
                // cout << "minimum_distance without size_chain = " << min_dist << endl;
                // cout << "maximun_distance without size_chain = " << max_dist << endl;
                // cout << "minimum_distance = " << minimum_distance << endl;
                // cout << "maximun_distance = " << maximun_distance << endl;

                is_complex = true;
            }
        }

        if (ppath.nreversed() > ppath.size() / 2) {
            ppath.flip();
        }

        if (is_complex) { // Case of complex found
            sum_path = 0;
        } else {
            for (size_t i = 1; i < size_node.size()-1; ++i) {
                sum_path += size_node[i];
            }
        }

        pretty_paths.push_back(ppath.print());
        size_t size_path = ppath.size();
        seq_net_paths.push_back(std::make_tuple(seq_net, minimum_distance, maximun_distance, size_path, sum_path));
    }

    vector<string> type_variants = calcul_pos_type_variant(seq_net_paths);
    return std::make_tuple(pretty_paths, type_variants);
}

// {chr : matrix(snarl, paths, start_pos, end_pos, type)}
std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>> loop_over_snarls_write(
        SnarlDistanceIndex& stree,
        vector<tuple<net_handle_t, string, size_t, size_t, bool>>& snarls,
        PackedGraph& pg, 
        const string& output_file,
        const string& output_snarl_not_analyse,
        size_t& children_threshold,
        size_t& path_length_threshold, 
        bool bool_return = true) {

    ofstream out_snarl(output_file);
    ofstream out_fail(output_snarl_not_analyse);
    
    out_snarl << "CHR\tPOS\tEND\tSNARL\tPATHS\tTYPE\tREF\n";
    out_fail << "SNARL\tREASON\n";
    
    std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>> snarl_paths;
    unordered_map<string, std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>> chr_snarl_matrix;
    size_t paths_number_analysis = 0;
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
            out_fail << snarl_id << "\ttoo_many_children = " << children[0] << " children" << "\n";
            continue;
        }
        
        vector<vector<net_handle_t>> paths = {{stree.get_bound(snarl, false, true)}};
        vector<vector<net_handle_t>> finished_paths;

        while (!paths.empty()) {
            auto path = paths.back();
            paths.pop_back();

            if (itr > path_length_threshold) {
                out_fail << snarl_id << "\titeration_calculation_out = " << children[0] << " children" << "\n";
                not_break = false;
                break;
            }
            follow_edges(stree, finished_paths, path, paths, pg);
            itr++;
        }

        if (not_break) {
            // pair<vector<string>, vector<string>>
            auto [pretty_paths, type_variants] = fill_pretty_paths(stree, pg, finished_paths);
            std::ostringstream pretty_paths_stream, type_variants_stream;

            // Convert pretty_paths (vector<string>) into a comma-separated string
            for (size_t i = 0; i < pretty_paths.size(); ++i) {
                if (i > 0) {pretty_paths_stream << ",";}
                pretty_paths_stream << pretty_paths[i];
            }
 
            for (size_t i = 0; i < type_variants.size(); ++i) {
                if (i > 0) {type_variants_stream << ",";}  // Add a comma and space between strings
                type_variants_stream << type_variants[i];
            }

            // snarl_id chromosome   start_position  end_position    paths    type
            string chr = std::get<1>(snarl_path_pos);
            size_t strat_pos = std::get<2>(snarl_path_pos);
            size_t end_pos = std::get<3>(snarl_path_pos);
            paths_number_analysis += pretty_paths.size();
            string str_reference = std::get<4>(snarl_path_pos) == true ? "0" : "1"; // 0 : on, 1 : out

            if (bool_return) {
                out_snarl << chr << "\t" << strat_pos << "\t" << end_pos
                    << "\t" << snarl_id << "\t" << pretty_paths_stream.str() 
                    << "\t" << type_variants_stream.str() << "\t" << str_reference << "\n";
            } else {
                // case new chr
                if (chr != save_chr && !save_chr.empty()) {
                    cout << "cleaning" << endl; 
                    chr_snarl_matrix[save_chr] = std::move(snarl_paths);
                    snarl_paths.clear();
                }
                save_chr = chr;
                snarl_paths.push_back(std::make_tuple(snarl_id, pretty_paths, strat_pos, end_pos, type_variants));
            }
        }
    }

    // last chr adding, but only if save_chr is not empty
    if (!save_chr.empty()) {
        chr_snarl_matrix[save_chr] = std::move(snarl_paths);
    }

    // Print the size of snarl_paths
    cout << "Number of paths : " << paths_number_analysis << endl;

    // Print chr_snarl_matrix
    for (const auto& chr_snarl : chr_snarl_matrix) {
        cout << "chr : " << chr_snarl.first << ", number of snarl : " << chr_snarl.second.size() << endl;
    }

    return {chr_snarl_matrix};
}
