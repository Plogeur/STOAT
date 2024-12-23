#include "list_snarl_paths.hpp"

using namespace std;
using namespace bdsg;

using Snarl = bdsg::SnarlDistanceIndex::NetHandle;

class Path {
private:
    std::vector<std::string> nodes;
    std::vector<char> orients;

public:
    // Constructor
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
};


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


vector<Snarl> save_snarls(SnarlDistanceIndex& stree, const Snarl& root) {
    vector<Snarl> snarls;

    function<void(const Snarl&)> save_snarl_tree_node = [&](const Snarl& net) {
        if (stree.is_snarl(net)) {
            snarls.push_back(net);
        }

        if (!stree.is_node(net) && !stree.is_sentinel(net)) {
            stree.for_each_child(net, save_snarl_tree_node);
        }
    };

    stree.for_each_child(root, save_snarl_tree_node);
    return snarls;
}

void parse_graph_tree(const string& pg_file, const string& dist_file, 
                      PackedGraph& pg, SnarlDistanceIndex& stree, Snarl& root) {
    pg.deserialize(pg_file);
    stree.deserialize(dist_file);
    root = stree.get_root();
}

pair<vector<string>, vector<vector<string>>> fill_pretty_paths(SnarlDistanceIndex& stree, 
                                                               PackedGraph& pg, 
                                                               const vector<vector<Snarl>>& finished_paths) {
    vector<string> pretty_paths;
    vector<vector<string>> length_net_paths;

    for (const auto& path : finished_paths) {
        string ppath;
        vector<string> length_net;

        for (const auto& net : path) {
            Snarl current_net = net;
            if (stree.is_sentinel(net)) {
                current_net = stree.get_node_from_sentinel(net);
            }

            if (stree.is_node(current_net) || stree.is_trivial_chain(current_net)) {
                // Add node handle or trivial chain
                if (stree.is_node(current_net)) {
                    length_net.push_back(to_string(stree.node_length(current_net)));
                } else {
                    auto stn_start = stree.get_bound(current_net, false, true);
                    auto node_start_id = stree.node_id(stn_start);
                    auto net_trivial_chain = pg.get_handle(node_start_id);
                    length_net.push_back(to_string(pg.get_length(net_trivial_chain)));
                }

            } else if (stree.is_chain(current_net)) {
                // Add chain bounds representation
                auto nodl = stree.get_bound(current_net, false, true);
                auto nodr = stree.get_bound(current_net, true, false);
                ppath += ">" + to_string(stree.node_id(nodl)) + ">*>" + to_string(stree.node_id(nodr));
                length_net.push_back("-1");
            }
        }

        pretty_paths.push_back(ppath);
        length_net_paths.push_back(length_net);
    }

    // Perform type variant calculation (placeholder function below)
    // length_net_paths = calcul_type_variant(length_net_paths);
    assert(length_net_paths.size() == pretty_paths.size());

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

void loop_over_snarls_write(SnarlDistanceIndex& stree, 
                            const vector<Snarl>& snarls, 
                            PackedGraph& pg, 
                            const string& output_file, 
                            const string& output_snarl_not_analyse, 
                            int time_threshold = 10) {
    write_header_output(output_file);

    for (const auto& snarl : snarls) {
        auto start_time = chrono::steady_clock::now();
        string snarl_id = to_string(stree.get_snarl_id(snarl)); // Placeholder

        vector<vector<Snarl>> paths = {{stree.get_bound(snarl, false, true)}};
        vector<vector<Snarl>> finished_paths;

        while (!paths.empty()) {
            auto path = paths.back();
            paths.pop_back();

            if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - start_time).count() > time_threshold) {
                // Log not analyzed snarls (if necessary)
                ofstream not_analyse_file(output_snarl_not_analyse, ios::app);
                not_analyse_file << snarl_id << "\ttime_calculation_out\n";
                break;
            }

            // follow_edges (Placeholder: Implement based on actual API requirements)
        }

        auto [pretty_paths, type_variants] = fill_pretty_paths(stree, pg, finished_paths);
        write_output(output_file, snarl_id, pretty_paths, type_variants);
    }
}
