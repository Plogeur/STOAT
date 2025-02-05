#include "list_snarl_paths.hpp"

using namespace std;
using Snarl = bdsg::SnarlDistanceIndex::NetHandle;

class Path {
private:
    std::vector<std::string> nodes;
    std::vector<char> orients;

public:
    Path() {}

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

string find_snarl_id(SnarlTree& stree, bdsg::SnarlTree::NetHandle snarl) {
    return stree.net_handle_as_string(snarl);
}

void follow_edges(SnarlTree& stree, vector<vector<bdsg::SnarlTree::NetHandle>>& finished_paths,
                  vector<bdsg::SnarlTree::NetHandle>& path,
                  vector<vector<bdsg::SnarlTree::NetHandle>>& paths, PackedGraph& pg) {
    // Implement edge following logic
}

pair<vector<tuple<string, vector<string>, string, string, string>>, int> loop_over_snarls_write(
        SnarlTree& stree, 
        vector<vector<bdsg::SnarlTree::NetHandle>> snarls,
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
        
        vector<vector<bdsg::SnarlTree::NetHandle>> paths = {{stree.get_bound(snarl, false, true)}};
        vector<vector<bdsg::SnarlTree::NetHandle>> finished_paths;
        
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
            
            if (bool_return) {
                snarl_paths.emplace_back(snarl_id, pretty_paths, join(type_variants, ","),
                                         snarl_path_pos[1], snarl_path_pos[2]);
            }
            
            paths_number_analysis += pretty_paths.size();
        }
    }
    
    return {snarl_paths, paths_number_analysis};
}

vector<Snarl> save_snarls(SnarlTree& stree, SnarlTree::NetHandle root,
                                   PackedGraph& pg, unordered_set<string>& ref_paths,
                                   PathPositionOverlay& ppo) {
    vector<Snarl> snarls;
    unordered_map<string, Position> snarls_pos;

    // Given a node handle (dist index) return a position on a reference path
    auto get_node_position = [&](SnarlTree::NetHandle node) -> Position {
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

    auto get_net_start_position = [&](bdsg::SnarlTree::NetHandle net) -> Position {
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

    function<void(bdsg::SnarlTree::NetHandle)> save_snarl_tree_node;
    save_snarl_tree_node = [&](bdsg::SnarlTree::NetHandle net) {
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
    bdsg::PackedGraph pg;
    bdsg::SnarlDistanceIndex stree;
    bdsg::PackedPositionOverlay pp_overlay;
    bdsg::net_handle_t root;

    // Constructor for GraphTree
    GraphTree(const std::string& pg_file, const std::string& dist_file) {
        // Initialize the attributes with deserialized data
        pg.deserialize(pg_file);
        stree.deserialize(dist_file);
        pp_overlay = bdsg::PackedPositionOverlay(pg);  // Initialize pp_overlay with the PackedGraph
        root = stree.get_root();  // Get the root from the SnarlDistanceIndex
    }

    // Getter for PackedGraph
    const bdsg::PackedGraph& get_pg() const {
        return pg;
    }

    // Getter for SnarlDistanceIndex
    const bdsg::SnarlDistanceIndex& get_stree() const {
        return stree;
    }

    // Getter for PackedPositionOverlay
    const bdsg::PackedPositionOverlay& get_pp_overlay() const {
        return pp_overlay;
    }

    // Getter for the root net handle
    bdsg::net_handle_t& get_root() const {
        return root;
    }
};

pair<vector<string>, vector<vector<string>>> fill_pretty_paths(
    SnarlTree& stree, PackedGraph& pg, 
    vector<vector<bdsg::SnarlTree::NetHandle>> finished_paths) {
    
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
                bdsg::SnarlTree::NetHandle nodl, nodr;
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

void print_help() {
    std::cout << "Usage: SnarlParser [options]\n\n"
              << "Options:\n"
              << "  -p, --pgfile <path>       Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -d, --distfile <path>          Path to the snarl file (.txt or .tsv)\n"
              << "  -o, --output <name>         Output name\n"
              << "  -h, --help                  Print this help message\n";
}


int main(int argc, char* argv[]) {

    bool show_help = false;
    std::string pg_file, dist_file;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-p" || arg == "--pgfile") && i + 1 < argc) {
            pg_file = argv[++i];
        } else if ((arg == "-d" || arg == "--distfile") && i + 1 < argc) {
            dist_file = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        } else {
            show_help = true;
        }
    }

    if (show_help || pg_file.empty() || dist_file.empty()) {
        print_help();
        return 0;
    }

    Tree = GraphTree(pg_file, dist_file);
    bdsg::PackedGraph pg = Tree.get_pg()
    bdsg::SnarlDistanceIndex stree = Tree.get_stree()
    bdsg::net_handle_t root = Tree.get_root()
    bdsg::PackedPositionOverlay pp_overlay = Tree.get_pp_overlay()

    snarls = save_snarls(stree, root, pg, reference, pp_overlay);
    SnarlDistanceIndex& stree, 
                            const vector<Snarl>& snarls, 
                            PackedGraph& pg, 
                            const string& output_file, 
                            const string& output_snarl_not_analyse, 
                            int time_threshold = 10
    loop_over_snarls_write(stree, snarls, pg, output_file)

    return EXIT_SUCCESS;
}