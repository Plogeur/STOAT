#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include "snarl_parser.hpp"     
#include "matrix.hpp"
#include "arg_parser.hpp"
#include "list_snarl_paths.hpp"

using namespace std;

void print_help() {
    std::cout << "Usage: SnarlParser [options]\n\n"
              << "Options:\n"
              << "  -v, --vcf_path <path>       Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -s, --snarl <path>          Path to the snarl file (.txt or .tsv)\n"
              << "  -p, --pg <path>             Path to the pg file (.pg)\n"
              << "  -d, --dist <path>           Path to the dist file (.dist)\n"
              << "  -b, --binary <path>         Path to the binary group file (.txt or .tsv)\n"
              << "  -q, --quantitative <path>   Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  -o, --output <name>         Output name\n"
              << "  -h, --help                  Print this help message\n";
}

int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, pg_path, dist_path, binary_path, quantitative_path, output_path;
    bool show_help = false;

    // Parse arguments manually
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-v" || arg == "--vcf_path") && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if ((arg == "-s" || arg == "--snarl") && i + 1 < argc) {
            snarl_path = argv[++i];
        } else if ((arg == "-p" || arg == "--pg") && i + 1 < argc) {
            pg_path = argv[++i];
        } else if ((arg == "-d" || arg == "--dist") && i + 1 < argc) {
            dist_path = argv[++i];
        } else if ((arg == "-b" || arg == "--binary") && i + 1 < argc) {
            binary_path = argv[++i];
        } else if ((arg == "-q" || arg == "--quantitative") && i + 1 < argc) {
            quantitative_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        }
    }

    std::filesystem::path output_dir = "output";
    std::filesystem::create_directory(output_dir);
    unordered_set<string> reference {"ref"};
    output_path = (output_dir / output_path).string();

    if (show_help || vcf_path.empty() ||
        binary_path.empty() == quantitative_path.empty()) {
        cerr << "vcf_path or phenotype are missing";
        print_help();
        return 0;
    }

    if (snarl_path.empty() && (pg_path.empty() || dist_path.empty())) {
        cerr << "snarl or pg and dist files are missing";
        print_help();
        return 0;
    }

    // Check format of the VCF file
    unordered_map<string, std::tuple<vector<string>, std::string, string, std::vector<string>>> snarl;

    if (!snarl_path.empty()){
        check_format_paths_snarl(snarl_path);
        snarl = parse_snarl_path(snarl_path);

    } else if (!pg_path.empty() && !dist_path.empty()) {
        auto [stree, pg, root, pp_overlay] = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, reference, *pp_overlay);
        string output_snarl_not_analyse = output_path + "/snarl_not_analyse.tsv";
        string output_file = output_path + "/snarl_analyse.tsv";
        int children_threshold = 50;
        snarl = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, true);
    
    } else {
        std::cerr << "Error: Either snarl file or pg and dist files must be provided\n";
        return EXIT_FAILURE;
    }

    std::vector<std::string> list_samples = parseHeader(vcf_path);    
    std::unordered_map<std::string, bool> binary;
    std::unordered_map<std::string, float> quantitative;

    if (!binary_path.empty()) {
        check_format_binary_phenotype(binary_path);
        binary = parse_binary_pheno(binary_path);
        check_match_samples(binary, list_samples);
    }

    if (!quantitative_path.empty()) {
        check_format_quantitative_phenotype(quantitative_path);
        quantitative = parse_quantitative_pheno(quantitative_path);
        check_match_samples(quantitative, list_samples);
    }

    // Initialize the SnarlProcessor with the VCF path
    SnarlParser vcf_object(vcf_path);
    auto start_1 = std::chrono::high_resolution_clock::now();
    vcf_object.fill_matrix();
    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time Matrix : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;
    auto start_2 = std::chrono::high_resolution_clock::now();

    // Process binary group file if provided
    if (!binary_path.empty()) {
        if (!output_path.empty()) {
            vcf_object.binary_table(snarl, binary, output_path);
        } else {
            vcf_object.binary_table(snarl, binary);
        }
    }

    // Process quantitative phenotype file if provided
    if (!quantitative_path.empty()) {
        if (!output_path.empty()) {
            vcf_object.quantitative_table(snarl, quantitative, output_path);
        } else {
            vcf_object.quantitative_table(snarl, quantitative);
        }
    }

    auto end_2 = std::chrono::high_resolution_clock::now();
    std::cout << "Time P-value : " << std::chrono::duration<double>(end_2 - start_2).count() << " s" << std::endl;
    std::cout << "Time Gwas analysis : " << std::chrono::duration<double>(end_2 - start_1).count() << " s" << std::endl;

    return EXIT_SUCCESS;
}

// ./stoat_cxx -p ../data/binary/pg.pg -d ../data/binary/pg.dist -v ../data/binary/binary.vcf.gz -b ../data/binary/phenotype.tsv -o output