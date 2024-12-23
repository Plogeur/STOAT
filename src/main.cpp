#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include "snarl_parser.hpp"     
#include "matrix.hpp"
#include "arg_parser.hpp"

void print_help() {
    std::cout << "Usage: SnarlParser [options]\n\n"
              << "Options:\n"
              << "  -v, --vcf_path <path>       Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -s, --snarl <path>          Path to the snarl file (.txt or .tsv)\n"
              << "  -b, --binary <path>         Path to the binary group file (.txt or .tsv)\n"
              << "  -q, --quantitative <path>   Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  -o, --output <name>         Output name\n"
              << "  -h, --help                  Print this help message\n";
}

int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, binary_path, quantitative_path, output_name;
    bool show_help = false;

    // Parse arguments manually
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-v" || arg == "--vcf_path") && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if ((arg == "-s" || arg == "--snarl") && i + 1 < argc) {
            snarl_path = argv[++i];
        } else if ((arg == "-b" || arg == "--binary") && i + 1 < argc) {
            binary_path = argv[++i];
        } else if ((arg == "-q" || arg == "--quantitative") && i + 1 < argc) {
            quantitative_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_name = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        }
    }

    std::filesystem::path output_dir = "output";
    std::filesystem::create_directory(output_dir);
    output_name = (output_dir / output_name).string();

    if (show_help || vcf_path.empty() || snarl_path.empty() || 
        binary_path.empty() == quantitative_path.empty()) {
        print_help();
        return 0;
    }

    // Check format of the VCF file
    check_format_vcf_file(vcf_path);

    // Check format of the snarl paths file
    check_format_paths_snarl(snarl_path);

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

    // Parse the snarl file
    auto snarl = parse_snarl_path(snarl_path);

    // Initialize the SnarlProcessor with the VCF path
    SnarlParser vcf_object(vcf_path);
    auto start_1 = std::chrono::high_resolution_clock::now();
    vcf_object.fill_matrix();
    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time Matrix: " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

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

    end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time P-value: " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

    return EXIT_SUCCESS;
}

// ./snarl_project --vcf_path ../test/small_vcf.vcf --snarl ../test/list_snarl_short.txt -b ../test/group.txt -o snarl_project
