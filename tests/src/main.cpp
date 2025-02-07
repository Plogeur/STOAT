#include "list_snarl_paths.hpp"

void print_help() {
    std::cout << "Usage: SnarlParser [options]\n\n"
              << "Options:\n"
              << "  -p, --pgfile <path>       Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -d, --distfile <path>     Path to the snarl file (.txt or .tsv)\n"
              << "  -o, --output <name>       Output path\n"
              << "  -h, --help                Print this help message\n";
}

int main(int argc, char* argv[]) {

    bool show_help = false;
    std::string pg_file, dist_file, output_path;

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

    if (show_help || pg_file.empty() || dist_file.empty() || output_path.empty()) {
        print_help();
        return 0;
    }

    unordered_set<string> reference {"ref"};
    std::cout << "Checking input files: " << pg_file << " " << dist_file << std::endl;
    auto [stree, pg, root, pp_overlay] = parse_graph_tree(pg_file, dist_file);
    std::cout << "Parsing successful!" << std::endl;
    vector<tuple<net_handle_t, string, size_t>> snarls = save_snarls(*stree, root, *pg, reference, *pp_overlay);
    std::cout << "saving snarls successful!" << std::endl;
    string output_snarl_not_analyse = output_path + "/snarl_not_analyse.tsv";
    string output_file = output_path + "/snarl_analyse.tsv";
    int children_threshold = 50;
    bool bool_return = true;
    loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, bool_return);
    std::cout << "Write snarls successful!" << std::endl;

    return EXIT_SUCCESS;
}

// ./list_snarl_paths -p ../../data/binary/pg.pg -d ../../data/binary/pg.dist -o ../binary_test
// terminate called after throwing an instance of 'std::invalid_argument'
//   what():  stoi
// Abandon (core dumped)