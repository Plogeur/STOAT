#include "snarl_data_t.hpp"

// using handlegraph::step_handle_t;
// using handlegraph::handle_t;
// using handlegraph::net_handle_t;

namespace stoat_vcf {

// Node_traversal_t
Node_traversal_t::Node_traversal_t(const size_t &id, const bool &rev)
        : node_id(id), is_reverse(rev) {}

// Convert Node_traversal_t to node + path representation [string]
std::string Node_traversal_t::to_string() const {
    return (is_reverse ? "<" : ">") + std::to_string(node_id);
}

// Getters for Node_traversal_t
size_t Node_traversal_t::get_node_id() const { return node_id; }
bool Node_traversal_t::get_is_reverse() const { return is_reverse; }

bool Node_traversal_t::operator==(const Node_traversal_t& other) const {
    return node_id == other.node_id && is_reverse == other.is_reverse;
}

// Edge_t
Edge_t::Edge_t(const Node_traversal_t &node_traversal_1, 
               const Node_traversal_t &node_traversal_2) :
    edge(std::make_pair(node_traversal_1, node_traversal_2)) {}

// Convert Edge_t to std::pair<size_t, size_t>
std::pair<size_t, size_t> Edge_t::print_pair_edge() const {
    return std::make_pair(edge.first.get_node_id(), edge.second.get_node_id());
}

// Convert Edge_t to std::string
std::string Edge_t::print_string_edge() const {
    return edge.first.to_string() + edge.second.to_string();
}

// Accessor to edge, useful for hashing and comparison
const std::pair<Node_traversal_t, Node_traversal_t>& Edge_t::get_edge() const {
    return edge;
}

bool Edge_t::operator==(const Edge_t &other) const {
    return edge == other.edge;
}

// add a node traversal to the path
void Path_traversal_t::add_node_traversal_t(const Node_traversal_t &node) {
    this->paths.push_back(node);
}

// convert Path_traversal_t to path representation
std::string Path_traversal_t::to_string() const {
    std::string result;
    for (const auto& node : paths) {
        result += node.to_string();
    }
    return result;
}

const std::vector<Node_traversal_t>& Path_traversal_t::get_paths() const { 
    return paths; 
};

std::string pairToString(const std::pair<size_t, size_t>& name) {
    std::ostringstream oss;
    oss << name.first << "_" << name.second;
    return oss.str();
}

std::pair<size_t, size_t> stringToPair(const std::string& str) {
    size_t underscorePos = str.find('_');
    if (underscorePos == std::string::npos) {
        throw std::invalid_argument("Input std::string does not contain an underscore separator");
    }

    std::string firstPart = str.substr(0, underscorePos);
    std::string secondPart = str.substr(underscorePos + 1);

    size_t first = std::stoul(firstPart);
    size_t second = std::stoul(secondPart);

    return {first, second};
}

std::string vectorPathToString(const std::vector<stoat_vcf::Path_traversal_t>& vec_paths) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec_paths.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec_paths[i].to_string();
    }
    return oss.str();
}

std::vector<stoat_vcf::Path_traversal_t> stringToVectorPath(std::string& input) {
    std::vector<stoat_vcf::Path_traversal_t> vec_paths;
    std::istringstream iss(input);
    std::string path_str;

    // Split by commas to get individual Path_traversal_t strings
    while (std::getline(iss, path_str, ',')) {
        Path_traversal_t path;
        size_t i = 0;

        while (i < path_str.size()) {
            // Parse node_id
            size_t node_id = 0;
            while (i < path_str.size() && std::isdigit(path_str[i])) {
                node_id = node_id * 10 + (path_str[i] - '0');
                ++i;
            }

            bool is_reverse = (path_str[i] == '<');
            ++i;

            Node_traversal_t node(node_id, is_reverse);
            path.add_node_traversal_t(node);
        }

        vec_paths.push_back(path);
    }

    return vec_paths;
}

// Add a snarl
Snarl_data_t::Snarl_data_t(bdsg::net_handle_t snarl_,
    std::vector<Path_traversal_t> snarl_paths_,
    const size_t start_positions_, const size_t end_positions_,
    std::vector<std::string> type_variants_) :
    snarl(snarl_),
    snarl_paths(std::move(snarl_paths_)),
    start_positions(start_positions_),
    end_positions(end_positions_),
    type_variants(std::move(type_variants_)) {}


Path::Path() {}

// Add a node with known orientation
void Path::addNode(const std::string& node, char orient) {
    nodes.push_back(node);
    orients.push_back(orient);
}

// Add a node handle and extract information using the std::string representation
bool Path::addNodeHandle(const handlegraph::net_handle_t& node_h, const bdsg::SnarlDistanceIndex& stree) {
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

// Get the std::string representation of the path
Path_traversal_t Path::print() const {
    Path_traversal_t out_path;
    for (size_t i = 0; i < nodes.size(); ++i) {
        size_t node_size_t;
        if (nodes[i] == "*") {
            node_size_t = 0; // Special case for "*""
        } else {
            node_size_t = std::stoi(nodes[i]);
        }
        Node_traversal_t node_traversal(
            node_size_t,
            orients[i] == '>' ? false : true // because is reverse is false for '>' and true for '<'
        );
        out_path.add_node_traversal_t(node_traversal);
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
// tuple<std::string, size_t, size_t, size_t>
// seq_net, minimum_distance, maximun_distance, size_path, sum_path
std::vector<std::string> calcul_pos_type_variant(const std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>>& list_length_paths) {
    std::vector<std::string> list_type_variant;

    for (const auto& tuple_info : list_length_paths) {
        size_t path_length = std::get<2>(tuple_info);
        size_t sum_path = std::get<3>(tuple_info);
        bool is_complex = std::get<4>(tuple_info);

        if (path_length >= 3) {
            if (is_complex) { // Case complex
                std::string complex = std::to_string(std::get<0>(tuple_info)) + "/" + std::to_string(std::get<1>(tuple_info));
                list_type_variant.push_back(complex);
            } else { // Case multiple nodes (ex : INS+SNP+...)
                list_type_variant.push_back(std::to_string(sum_path));
            }

        } else if (path_length == 2) { // case Deletion
            list_type_variant.push_back("0");
        } else { // Case path_lengths is empty or == 1
            std::cerr << "path_lengths is empty" << std::endl;
        }
    }
    return list_type_variant;
}

std::pair<size_t, size_t> find_snarl_id(const bdsg::SnarlDistanceIndex& stree, const handlegraph::net_handle_t& snarl) {
    
    // Get start and end boundary nodes for the snarl
    auto sstart = stree.get_bound(snarl, false, true);  // False for the left boundary
    auto send = stree.get_bound(snarl, true, true);     // True for the right boundary

    // Convert the sentinels into nodes
    auto start_node = stree.get_node_from_sentinel(sstart);
    auto end_node = stree.get_node_from_sentinel(send);

    // Get the node IDs from bdsg::SnarlDistanceIndex
    // handlegraph::nid_t
    auto start_node_id = stree.node_id(start_node);
    auto end_node_id = stree.node_id(end_node);

    // Convert to size_t
    size_t start_node_id_size_t = static_cast<size_t>(start_node_id);
    size_t end_node_id_size_t = static_cast<size_t>(end_node_id);

    // Construct the snarl ID
    std::pair<size_t, size_t> snarl_id(end_node_id_size_t, start_node_id_size_t);

    return snarl_id;  // Return the generated snarl ID as a std::string
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

    //bdsg::PackedPositionOverlay takes a pointer to pg
    auto pp_overlay = std::make_unique<bdsg::PackedPositionOverlay>(pg.get());

    // Get root of snarl tree
    handlegraph::net_handle_t root = stree->get_root();

    return std::make_tuple(std::move(stree), std::move(pg), root, std::move(pp_overlay));
}

void follow_edges(bdsg::SnarlDistanceIndex& stree,
                std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths,
                const std::vector<handlegraph::net_handle_t>& path,
                std::vector<std::vector<handlegraph::net_handle_t>>& paths,
                bdsg::PackedGraph& pg,
                const bool& cycle) {

    auto add_to_path = [&](const handlegraph::net_handle_t& next_child) {

        // If this is the bound of the snarl then we're done && next_child is different that the first node
        if (stree.is_sentinel(next_child)) {
            size_t next_child_node_id = stree.node_id(stree.get_node_from_sentinel(next_child));
            size_t first_element_path_node_id = stree.node_id(stree.get_node_from_sentinel(path[0]));
            if (next_child_node_id != first_element_path_node_id) {
                finished_paths.emplace_back(path);
                finished_paths.back().push_back(next_child);
            }

        } else {

            if (cycle) { // Case where we find a loop
                return false;
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

std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>> save_snarls(
                                bdsg::SnarlDistanceIndex& stree, 
                                handlegraph::net_handle_t& root,
                                bdsg::PackedGraph& pg, 
                                std::unordered_set<std::string>& ref_chr,
                                bdsg::PackedPositionOverlay& ppo) {

    std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>> snarls;
    unordered_map<std::string, std::tuple<std::string, size_t, size_t>> snarls_pos;
    size_t save_end_pos_ref = 0;

    // Given a node handle (dist index), return a position if on chr reference path
    auto get_node_position = [&](handlegraph::net_handle_t node) -> std::tuple<std::string, size_t, size_t> { // node : handlegraph::net_handle_t
        handlegraph::handle_t node_h = stree.get_handle(node, &pg);

        // path_name, position
        std::tuple<std::string, size_t, size_t> ret_pos;

        auto step_callback = [&](const handlegraph::step_handle_t& step_handle) {
            handlegraph::path_handle_t path_handle = pg.get_path_handle_of_step(step_handle);
            std::string chr_path = pg.get_path_name(path_handle);

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

    auto get_net_start_position = [&](handlegraph::net_handle_t net) -> std::tuple<std::string, size_t, size_t> {

        if (stree.is_node(net)) {
            return get_node_position(net);
        }

        handlegraph::net_handle_t bnode1 = stree.get_bound(net, true, false);
        std::tuple<std::string, size_t, size_t> bnode1_p = get_node_position(bnode1);

        handlegraph::net_handle_t bnode2 = stree.get_bound(net, false, false); // verify false true ?
        std::tuple<std::string, size_t, size_t> bnode2_p = get_node_position(bnode2);

        // Check if the std::string part of the pair is empty
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

        // tuple<std::string, size_t, size_t> snarl_start_end;
        return make_tuple(std::get<0>(bnode1_p), start, end);
    };

    function<void(handlegraph::net_handle_t)> save_snarl_tree_node;
    save_snarl_tree_node = [&](handlegraph::net_handle_t net) {

        std::tuple<std::string, size_t, size_t> snarl_pos = get_net_start_position(net);
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
            // handlegraph::net_handle_t snarl, chr_ref, pos, is_on_ref_bool
            snarls.push_back(std::make_tuple(net, std::get<0>(snarl_pos), std::get<1>(snarl_pos), std::get<2>(snarl_pos), bool_ref));
        }

        // explore children
        if (!stree.is_node(net) && !stree.is_sentinel(net)) {
            stree.for_each_child(net, save_snarl_tree_node);
        }
    };

    stree.for_each_child(root, save_snarl_tree_node);
    cout << "Number of snarls : " << snarls.size() << std::endl;
    return snarls;
}

std::tuple<std::vector<stoat_vcf::Path_traversal_t>, std::vector<std::string>> fill_pretty_paths(
    bdsg::SnarlDistanceIndex& stree, 
    bdsg::PackedGraph& pg, 
    std::vector<std::vector<handlegraph::net_handle_t>>& finished_paths) {
    
    // list of paths
    std::vector<stoat_vcf::Path_traversal_t> pretty_paths;

    // seq_net, minimum_distance, maximun_distance, size_path, sum_path
    // Used to calculate the type of variant
    std::vector<std::tuple<size_t, size_t, size_t, size_t, bool>> seq_net_paths;

    for (const auto& path : finished_paths) {
        Path ppath;
        bool is_complex = false;
        size_t sum_path = 0;
        size_t minimum_distance=0;
        size_t maximun_distance=0;
        std::vector<size_t> size_node;
        size_node.resize(path.size(), 0);

        for (int i=0; i<path.size(); i++) {
            handlegraph::net_handle_t net = path[i];

            if (stree.is_sentinel(net)) {
                net = stree.get_node_from_sentinel(net);
            }

            // Node case
            if (stree.is_node(net)) {
                bool rev = ppath.addNodeHandle(net, stree);
                handlegraph::nid_t node_start_id = stree.node_id(net);
                handlegraph::handle_t node_handle = pg.get_handle(node_start_id);
                size_node[i] = pg.get_length(node_handle);
            }

            // Trivial chain case
            else if (stree.is_trivial_chain(net)) {
                bool rev = ppath.addNodeHandle(net, stree);
                auto stn_start = stree.starts_at_start(net) ? stree.get_bound(net, false, true) : stree.get_bound(net, true, true);
                handlegraph::nid_t node_start_id = stree.node_id(stn_start);
                handlegraph::handle_t net_trivial_chain = pg.get_handle(node_start_id);
                size_node[i] = pg.get_length(net_trivial_chain);
            }

            // Chain case aka complex
            else if (stree.is_chain(net)) {
                handlegraph::net_handle_t nodl, nodr;
                if (stree.starts_at_start(net)) {
                    nodl = stree.get_bound(net, false, true);
                    nodr = stree.get_bound(net, true, false);
                } else {
                    nodl = stree.get_bound(net, true, true);
                    nodr = stree.get_bound(net, false, false);
                }

                ppath.addNodeHandle(nodl, stree);

                // TODO? test chain : handlegraph::net_handle_t is composed of 2 element && if both element is_node == true ?
                // idk ask to jean
                bool chain_2node = true;
                int child_count = 0;
                size_t sum_node = 0;

                stree.for_each_child(net, [&](const handlegraph::net_handle_t& child) {
                    ++child_count;
                    if (!stree.is_node(child)) {
                        chain_2node = false;
                        return false; // stop early
                    } else {
                        sum_node += pg.get_length(pg.get_handle(stree.node_id(child)));
                    }
                    return true;
                });
                
                if (!(chain_2node && child_count == 2)) {
                    ppath.addNode("*", '>');
                    is_complex = true;
                } else {
                    size_node[i] = sum_node;
                }
                ppath.addNodeHandle(nodr, stree);

                // Get the size of the chain and return the distance (minimum and maximum)
                size_t complex_start_id = stree.node_id(nodl);
                handlegraph::handle_t handle_start = pg.get_handle(complex_start_id);
                size_t size_start_node = pg.get_length(handle_start);
                bool revl = stree.ends_at_start(nodl);

                size_t complex_end_id = stree.node_id(nodr);
                handlegraph::handle_t handle_end = pg.get_handle(complex_end_id);
                size_t size_end_node = pg.get_length(handle_end);
                bool revr = stree.ends_at_start(nodr);

                size_t size_chain = size_start_node + size_end_node;
                // TODO: I think this can use minimum_length() and maximum_length(), just to be simpler
                // matis ans : yes for minimum_length() but maximum_length() do not exist 
                size_t min_dist = stree.minimum_distance(complex_start_id, revl, size_start_node, complex_end_id, revr, 0);
                size_t max_dist = stree.maximum_distance(complex_start_id, revl, size_start_node, complex_end_id, revr, 0);

                // Fail case 
                assert(max_dist != static_cast<size_t>(INT_MAX) && "Overflow max distance");
                assert(min_dist != static_cast<size_t>(INT_MAX) && "Overflow min distance");

                minimum_distance += size_chain + min_dist;
                maximun_distance += size_chain + max_dist;
            }
        }

        if (ppath.nreversed() > ppath.size() / 2) {
            ppath.flip();
        }

        if (is_complex) { // Case of complex found
            for (size_t i = 1; i < size_node.size()-1; ++i) {
                maximun_distance += size_node[i];
                minimum_distance += size_node[i];
            }
        } else {
            for (size_t i = 1; i < size_node.size()-1; ++i) {
                sum_path += size_node[i];
            }
        }

        pretty_paths.push_back(ppath.print());
        size_t size_path = ppath.size();
        seq_net_paths.push_back(std::make_tuple(minimum_distance, maximun_distance, size_path, sum_path, is_complex));
    }

    // TODO : change sum_path to use boundary to compute it / remove size_node_2 
    // Matis ans : I remove size_node_2 BUT i don't know how to use boundary to compute it.
    std::vector<std::string> type_variants = calcul_pos_type_variant(seq_net_paths);
    return std::make_tuple(pretty_paths, type_variants);
}

// {chr : matrix(snarl, paths, start_pos, end_pos, type)}
std::unordered_map<std::string, std::vector<Snarl_data_t>> loop_over_snarls_write(
        bdsg::SnarlDistanceIndex& stree,
        std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>>& snarls,
        bdsg::PackedGraph& pg, 
        const std::string& output_file,
        const std::string& output_snarl_not_analyse,
        const size_t& children_threshold,
        const size_t& path_length_threshold, 
        const size_t& cycle_threshold,
        bool bool_return = true) {

    ofstream out_snarl(output_file);
    if (bool_return) {
        out_snarl << "CHR\tSTART_POS\tEND_POS\tSNARL\tPATHS\tTYPE\tREF\n";
    }

    ofstream out_fail(output_snarl_not_analyse);
    out_fail << "SNARL\tREASON\n";
        
    std::vector<Snarl_data_t> snarl_paths;
    snarl_paths.reserve(snarls.size()); // Reserve snarls size to avoid reallocations

    unordered_map<std::string, std::vector<Snarl_data_t>> chr_snarl_matrix;
    size_t paths_number_analysis = 0;
    std::string save_chr = "";

    // TODO: I think this should just be a size_t child_count, it doesn't look like it ever uses anything except the first value
    // Matis ans : i think size_t variable isn't in the scope of the lambda, so i use a vector with one element
    // but i agree that it is not the best way to do it, i will change if it's work with size_t
    std::vector<size_t> children = {0};
    auto count_children = [&](handlegraph::net_handle_t net) {
        children[0] += 1;
        return true;
    };

    for (const auto& snarl_path_pos : snarls) {
        handlegraph::net_handle_t snarl = std::get<0>(snarl_path_pos);
        size_t itr = 0;
        std::string snarl_id_str = pairToString(find_snarl_id(stree, snarl));
        bool not_break = true;
        children = {0}; // re-initialise the children vec
        
        stree.for_each_child(snarl, count_children);
        if (children[0] > children_threshold) {
            out_fail << snarl_id_str << "\ttoo_many_children = " << children[0] << " children" << "\n";
            continue;
        }
        
        // Find all paths going through the netgraph of the snarl
        std::vector<std::vector<handlegraph::net_handle_t>> paths = {{stree.get_bound(snarl, false, true)}};
        std::vector<std::vector<handlegraph::net_handle_t>> finished_paths;

        while (!paths.empty()) {
            //TODO: I think path should be a reference so it doesn't get copied
            //Matis ans: No it elements must remain in copy because it will be modified later (i test the & and it breaks the code : 0 paths found)
            //Change to move instead
            std::vector<handlegraph::net_handle_t> path = paths.back();
            std::unordered_map<handlegraph::net_handle_t, size_t> dict_path_occ;
            bool cycle = false;

            for (const auto& net : path) {
                dict_path_occ[net]++;
                if (dict_path_occ[net] > cycle_threshold+1) {
                    cycle = true;
                    break;
                }
            }

            paths.pop_back();

            if (itr > path_length_threshold) {
                out_fail << snarl_id_str << "\titeration_calculation_out = " << children[0] << " children" << "\n";
                not_break = false;
                break;
            }
            follow_edges(stree, finished_paths, path, paths, pg, cycle);
            itr++;
        }

        if (not_break) {
            // pair<std::vector<std::string>, std::vector<std::string>>
            auto [pretty_paths, type_variants] = fill_pretty_paths(stree, pg, finished_paths);

            std::string chr = std::get<1>(snarl_path_pos);
            size_t pretty_paths_size = pretty_paths.size();

            // skip this snarl with no chr ref in it OR with snarl less than 2 paths
            if (chr.empty() || pretty_paths_size < 2) {
                continue;
            }

            size_t strat_pos = std::get<2>(snarl_path_pos);
            size_t end_pos = std::get<3>(snarl_path_pos);
            paths_number_analysis += pretty_paths_size;
            std::string str_reference = std::get<4>(snarl_path_pos) == true ? "1" : "0"; // 1 : on reference, 0 : out reference

            if (bool_return) {
                out_snarl << chr << "\t" << strat_pos << "\t" << end_pos
                    << "\t" << handlegraph::as_integer(snarl) << "\t" << vectorPathToString(pretty_paths)
                    << "\t" << vectorToString(type_variants) << "\t" << str_reference << "\n";
            } else {
                // case new chr
                if (chr != save_chr && !save_chr.empty()) {
                    chr_snarl_matrix[save_chr] = std::move(snarl_paths);
                    snarl_paths.clear();
                }
                save_chr = chr;
                Snarl_data_t snarl_path(snarl, pretty_paths, strat_pos, end_pos, type_variants);
                snarl_paths.push_back(snarl_path);
            }
        }
    }

    // last chr adding, but only if save_chr is not empty
    if (!save_chr.empty()) {
        chr_snarl_matrix[save_chr] = std::move(snarl_paths);
    }

    // Print the size of snarl_paths
    cout << "Number of paths : " << paths_number_analysis << std::endl;

    // Print chr_snarl_matrix
    for (const auto& chr_snarl : chr_snarl_matrix) {
        cout << "chr : " << chr_snarl.first << ", number of snarl : " << chr_snarl.second.size() << std::endl;
    }

    return {chr_snarl_matrix};
}

} //end stoat_vcf namespace

// vg find -x ../snarl_data/fly.gbz -r 5176878:5176884 -c 10 | vg view -dp - | dot -Tsvg -o ../snarl_data/subgraph.svg
