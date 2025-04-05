#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <Eigen/Dense>
#include <cstdlib>

#include "snarl_parser.hpp"     
#include "matrix.hpp"
#include "arg_parser.hpp"
#include "list_snarl_paths.hpp"
#include "gaf_creator.hpp"

using namespace std;

void print_help() {
    std::cout << "Usage: SnarlParser [options]\n\n"
              << "Options:\n"
              << "  -v, --vcf_path <path>       Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -s, --snarl <path>          Path to the snarl file (.txt or .tsv)\n"
              << "  -p, --pg <path>             Path to the pg file (.pg)\n"
              << "  -d, --dist <path>           Path to the dist file (.dist)\n"
              << "  -r, --chr_ref <path>        Path to the chromosome reference file (.txt)\n"
              << "  -c, --children <int>        Max number of children for a snarl in the snarl decomposition process (default = 50)\n"
              << "  -b, --binary <path>         Path to the binary group file (.txt or .tsv)\n"
              << "  -g, --gaf                   Make a GAF file from the GWAS analysis\n"
              << "  -q, --quantitative <path>   Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  --covariate <path>          Path to the covariate file (.txt or .tsv)\n"
              << "  -e, --eqtl <path>           Path to the Expression Quantitative Trait Loci file (.txt or .tsv)\n"
              << "  --make-bed                  Create a plink format files (.bed, .bim, bed)\n"
              << "  --plot                      Create Manhatthan plot and QQ plot\n"
              << "  --maf                       Add a maf (Maximum allele frequency) thresold (defauld : 0.99)\n"
              << "  -o, --output <name>         Output dir name\n"
              << "  -t, --thread <int>          Number of threads\n"
              << "  -h, --help                  Print this help message\n";
}

int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, pg_path, dist_path, 
        chromosome_path, binary_path, quantitative_path, 
        eqtl_path, covariate_path, output_dir;

    size_t threads=1;
    size_t phenotype=0;
    size_t children_threshold = 50;
    double maf = 0.99;
    bool gaf = false;
    bool only_snarl_parsing = false;
    bool show_help = false;
    bool make_bed = false;
    bool make_plot = false;

    // Parse arguments manually
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-v" || arg == "--vcf_path") && i + 1 < argc) {
            vcf_path = argv[++i];
            check_file(vcf_path);
        } else if ((arg == "-s" || arg == "--snarl") && i + 1 < argc) {
            snarl_path = argv[++i];
            check_file(snarl_path);
        } else if ((arg == "-p" || arg == "--pg") && i + 1 < argc) {
            pg_path = argv[++i];
            check_file(pg_path);
        } else if ((arg == "-d" || arg == "--dist") && i + 1 < argc) {
            dist_path = argv[++i];
            check_file(dist_path);
        } else if ((arg == "-r" || arg == "--chr_ref") && i + 1 < argc) {
            chromosome_path = argv[++i];
            check_file(chromosome_path);
        } else if ((arg == "--make-bed") && i + 1 < argc) {
            make_bed=true;
        } else if ((arg == "--plot") && i + 1 < argc) {
            make_plot=true;
        } else if ((arg == "-c" || arg == "--children") && i + 1 < argc) {
            children_threshold = std::stoi(argv[++i]);
            if (children_threshold < 2) {
                std::cerr << "Error: Number of children must be a positive integer > 1\n";
                return EXIT_FAILURE;
            }
        } else if ((arg == "-b" || arg == "--binary") && i + 1 < argc) {
            binary_path = argv[++i];
            phenotype ++;
            check_file(binary_path);
        } else if ((arg == "-g" || arg == "--gaf") && i + 1 < argc) {
            gaf=true;
        } else if ((arg == "--covariate") && i + 1 < argc) {
            covariate_path = argv[++i];
            check_file(covariate_path);
        } else if ((arg == "-q" || arg == "--quantitative") && i + 1 < argc) {
            quantitative_path = argv[++i];
            phenotype ++;
            check_file(quantitative_path);
        } else if ((arg == "-e" || arg == "--eqtl") && i + 1 < argc) {
            eqtl_path = argv[++i];
            phenotype ++;
            check_file(eqtl_path);
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            // convert str to int and verify that it is a positive number
            threads = std::stoi(argv[++i]);
            if (threads < 1) {
                std::cerr << "Error: Number of threads must be a positive integer\n";
                return EXIT_FAILURE;
            }
        } else if ((arg == "--maf") && i + 1 < argc) {
            // convert str to int and verify that it is a positive number
            maf = std::stoi(argv[++i]);
            if (maf < 0 || maf > 1) {
                std::cerr << "Error: maf threshold must be a between 0 and 1\n";
                return EXIT_FAILURE;
            }
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            output_dir = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            show_help = true;
        }
    }

    if (output_dir.empty()) {
        output_dir = "../output";
    }

    auto start_1 = std::chrono::high_resolution_clock::now();
    std::filesystem::create_directory(output_dir);
    std::unordered_set<std::string> ref_chr = (!chromosome_path.empty()) ? parse_chromosome_reference(chromosome_path) : std::unordered_set<std::string>{"ref"};

    if (show_help) {
        print_help();
        return EXIT_FAILURE;    
    }

    // Enforce valid argument combinations
    if ((!snarl_path.empty() || (!pg_path.empty() && !dist_path.empty())) && !vcf_path.empty() && phenotype == 1) {
        // Case 1: snarl_path + vcf_path + phenotype
        // Case 2: pg_path + dist_path + vcf_path + phenotype
    } else if (!pg_path.empty() && !dist_path.empty() && vcf_path.empty() && snarl_path.empty() && phenotype == 0) {
        // Case 3: Only pg_path + dist_path
        only_snarl_parsing = true;
    } else if (((!pg_path.empty() && !dist_path.empty()) || (!snarl_path.empty())) && !vcf_path.empty() && make_bed == true) {
        // Case 4: Only pg_path + dist_path + vcf_path + make_bed activated
        // Case 5: snarl_path + vcf_path + --make-bed
    } else {
        std::cerr << "Invalid argument combination provided.\n";
        std::cerr << "There are 5 ways to lauch stoat : " << endl;
        std::cerr << "Case 1: snarl_path + vcf_path + phenotype (+ optional file)" << endl;
        std::cerr << "Case 2: pg_path + dist_path + vcf_path + phenotype (+ optional file)" << endl;
        std::cerr << "Case 3: pg_path + dist_path" << endl;
        std::cerr << "Case 4: pg_path + dist_path + vcf_path + --make-bed" << endl;
        std::cerr << "Case 5: snarl_path + vcf_path + --make-bed" << endl;

        print_help();
        return EXIT_FAILURE;
    }

    if ((gaf == true && binary_path.empty()) || (gaf == true && pg_path.empty())) {
        cerr << "GAF file can be generated only with binary phenotype AND with the pg graph";
        print_help();
        return EXIT_FAILURE;
    }

    // Check phenotypes
    std::vector<std::string> list_samples;
    htsFile* ptr_vcf;
    bcf_hdr_t* hdr;
    bcf1_t* rec;

    if (!only_snarl_parsing) {
        std::tie(list_samples, ptr_vcf, hdr, rec) = parseHeader(vcf_path);    
    }

    std::unordered_map<std::string, bool> binary;
    std::unordered_map<std::string, double> quantitative;
    //std::vector<QTLRecord> eqtl;

    std::unordered_map<std::string, std::vector<double>> covariate;
    if (!covariate_path.empty()) {
        check_format_covariate(covariate_path);
        covariate = parseCovariate(covariate_path);
    }

    if (!binary_path.empty()) {
        check_format_binary_phenotype(binary_path);
        binary = parse_binary_pheno(binary_path);
        check_match_samples(binary, list_samples);
        if (!covariate_path.empty()) {
            check_phenotype_covariate(binary, covariate);
        }

    } else if (!quantitative_path.empty()) {
        check_format_quantitative_phenotype(quantitative_path);
        quantitative = parse_quantitative_pheno(quantitative_path);
        check_match_samples(quantitative, list_samples);
        if (!covariate_path.empty()) {
            check_phenotype_covariate(quantitative, covariate);
        }

    } else if (!eqtl_path.empty()) {
        //eqtl = parseQTLFile(eqtl_path);
        //check_match_samples_eqtl(eqtl, list_samples);
    }

    // scope declaration
    // chr : <snarl, paths, pos, type>
    std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> snarls_chr;
    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<bdsg::PackedGraph> pg;
    handlegraph::net_handle_t root;
    std::unique_ptr<bdsg::PackedPositionOverlay> pp_overlay;
    size_t total_snarl;

    if (!snarl_path.empty()){
        std::tie(snarls_chr, total_snarl) = parse_snarl_path(snarl_path);
    } else {
        std::cout << "Start snarl analysis... " << std::endl;
        auto start_0 = std::chrono::high_resolution_clock::now();
        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);

        string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
        string output_file = output_dir + "/snarl_analyse.tsv";

        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, only_snarl_parsing);
        auto end_0 = std::chrono::high_resolution_clock::now();
        std::cout << "Snarl analysis : " << std::chrono::duration<double>(end_0 - start_0).count() << " s" << std::endl;
        if (only_snarl_parsing) {
            return EXIT_SUCCESS;
        }

        // Clean up unique_ptr except pg
        stree.reset();
        pp_overlay.reset();
    }

    if (make_bed) {
        std::vector<std::pair<std::string, int>> pheno;
        
        for (const auto& sample : list_samples) {
            pheno.push_back({sample, -9}); // initilize all phenotypes to -9
        }

        const std::string output_fam = output_dir + "/genotype.fam";
        create_fam(pheno, output_fam);
        // chromosome_chuck_make_bed(ptr_vcf, hdr, rec, list_samples, snarls_chr, output_dir);

        auto end_1 = std::chrono::high_resolution_clock::now();
        std::cout << "Time genotype plink files creations : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;
        return EXIT_SUCCESS;

    } else if (!binary_path.empty()) {

        string output_binary = output_dir + "/binary_analysis.tsv";
        chromosome_chuck_binary(ptr_vcf, hdr, rec, list_samples, snarls_chr, binary, covariate, maf, total_snarl, output_binary);

        string output_significative = output_dir + "/top_variant_binary.tsv";
        std::string python_cmd = "python3 ../src/adjusted_pvalue.py "
        " --pvalue " + output_binary + 
        " --binary";
        system(python_cmd.c_str());

        if (gaf) {
            string output_gaf = output_dir + "/binary_analysis.gaf";
            gaf_creation(output_binary, snarls_chr, *pg, output_gaf);
        }

        if (make_plot) {
            string output_manh = output_dir + "/manhattan_plot_binary.png";
            string output_qq = output_dir + "/qq_plot_binary.png";
            std::string python_cmd = "python3 ../src/make_plots.py "
            " --pvalue " + output_binary + 
            " --qq " + output_qq + 
            " --manh " + output_manh +
            " --binary";
            system(python_cmd.c_str());
        }

    } else if (!quantitative_path.empty()) {

        string output_quantitive = output_dir + "/quantitative_analysis.tsv";
        chromosome_chuck_quantitative(ptr_vcf, hdr, rec, list_samples, snarls_chr, quantitative, covariate, maf, total_snarl, output_quantitive);
        
        string output_significative = output_dir + "/top_variant_quantitative.tsv";
        std::string python_cmd = "python3 ../src/adjusted_pvalue.py "
        " --pvalue " + output_quantitive + 
        " --quantitative";
        system(python_cmd.c_str());

        if (make_plot) {
            string output_manh = output_dir + "/manhattan_plot_quantitative.png";
            string output_qq = output_dir + "/qq_plot_quantitative.png";

            std::string python_cmd = "python3 ../src/make_plots.py "
            " --pvalue " + output_quantitive + 
            " --qq " + output_qq + 
            " --manh " + output_manh +
            " --quantitative";
            system(python_cmd.c_str());
        }

    } else if (!eqtl_path.empty()) {
        string eqtl_output = output_dir + "/eqtl_gwas.tsv";
        std::ofstream outf(eqtl_output, std::ios::binary);
        std::string headers = "CHR\tPOS\tSNARL\tTYPE\tSE\tBETA\tP\n";
        outf.write(headers.c_str(), headers.size());

        // chromosome_chuck_eqtl(ptr_vcf, hdr, rec, list_samples, snarls_chr, eqtl, outf);
    }
    
    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time Gwas analysis : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

    return EXIT_SUCCESS;
}

// BINARY
// ./stoat_cxx -p ../data/binary/pg.pg -d ../data/binary/pg.dist -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv

// QUANTITATIVE
// ./stoat_cxx -p ../data/quantitative/pg.pg -d ../data/quantitative/pg.dist -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv
