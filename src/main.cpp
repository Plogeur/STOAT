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
              << "  --pheno <path>              Path to the binary group file (.txt or .tsv)\n"
              << "  --sex <path>                Make a plink format base on the snarl paths list and vcf file\n"
              << "  -o, --output <prefix>       Output prefix name\n"
              << "  -h, --help                  Print this help message\n";
}

int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, pheno_path, sex_path, prefix_name;
    bool show_help = false;

    // Parse arguments manually
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-v" || arg == "--vcf_path") && i + 1 < argc) {
            vcf_path = argv[++i];
        } else if ((arg == "-s" || arg == "--snarl") && i + 1 < argc) {
            snarl_path = argv[++i];
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            prefix_name = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        }
    }

    std::filesystem::path output_dir = "../output";
    std::filesystem::create_directory(output_dir);
    if (prefix_name.empty()) {prefix_name = "genotype";}
    const std::string output_name = (output_dir / prefix_name).string();

    if (show_help || vcf_path.empty() || snarl_path.empty()) {
        print_help();
        return 0;
    }

    // Check format of the VCF file
    check_format_vcf_file(vcf_path);

    // Check format of the snarl paths file
    check_format_paths_snarl(snarl_path);

    std::vector<std::string> list_samples = parseHeader(vcf_path);
 
    // Parse the snarl file
    auto snarl = parse_snarl_path(snarl_path);

    // Initialize the SnarlProcessor with the VCF path
    SnarlParser vcf_object(vcf_path);
    auto start_1 = std::chrono::high_resolution_clock::now();
    vcf_object.fill_matrix();
    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time matrix: " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

    // create the bim/bed plink format
    const std::string output_bim = output_name + ".bim";
    const std::string output_bed = output_name + ".bed";
    const std::string output_fam = output_name + ".fam";
    
    create_fam(list_samples, output_fam);
    vcf_object.create_bim_bed(snarl, output_bim, output_bed);

    end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time plink files creations : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

    return EXIT_SUCCESS;
}

// ./slink --vcf_path ../data/binary/merged_output.vcf --snarl ../data/binary/snarl_paths.tsv --pheno ../data/binary/phenotype.tsv -o genotype
