#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
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
              << "  -r, --chr_ref <path>        Path to the chromosome reference file (.tsv)\n"
              << "  -c, --children <int>        Max number of children for a snarl in the snarl decomposition process\n"
              << "  -b, --binary <path>         Path to the binary group file (.txt or .tsv)\n"
              << "  -g, --gaf                   Make GAF file from the GWAS analysis\n"
              << "  -q, --quantitative <path>   Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  -e, --eqtl <path>           Path to the Expression Quantitative Trait Loci file (.txt or .tsv)\n"
              << "  -o, --output <name>         Output dir name\n"
              << "  -t, --thread <int>          Number of threads\n"
              << "  -h, --help                  Print this help message\n";
}

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
                        const std::vector<std::string> &list_samples,
                        unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
                        const unordered_map<string, double>& pheno, std::ofstream& outf) {

    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << "GWAS analysis for chromosome : " << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.quantitative_table(snarl, pheno, chr, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
                        const std::vector<std::string> &list_samples, 
                        unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
                        const unordered_map<string, bool>& pheno, std::ofstream& outf) {

    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << "GWAS analysis for chromosome : " << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.binary_table(snarl, pheno, chr, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, pg_path, dist_path, chromosome_path, binary_path, quantitative_path, eqtl_path, output_dir;
    size_t threads=1;
    size_t phenotype=0;
    bool gaf, show_help= false;
    size_t children_threshold = 50;

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
        } else if ((arg == "-c" || arg == "--children") && i + 1 < argc) {
            children_threshold = std::stoi(argv[++i]);
            if (children_threshold < 1) {
                std::cerr << "Error: Number of children must be a positive integer\n";
                return EXIT_FAILURE;
            }
        } else if ((arg == "-b" || arg == "--binary") && i + 1 < argc) {
            binary_path = argv[++i];
            phenotype ++;
            check_file(binary_path);
        } else if ((arg == "-g" || arg == "--gaf") && i + 1 < argc) {
            gaf=true;
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

    if (vcf_path.empty()) {
        cerr << "vcf_path are missing";
        print_help();
        return EXIT_FAILURE;
    }

    if (phenotype == 0) {
        cerr << "phenotype are missing (use -b or -q or -eqtl)";
        print_help();
        return EXIT_FAILURE;
    }

    if (snarl_path.empty() && (pg_path.empty() || dist_path.empty())) {
        cerr << "snarl paths file OR pg & dist files are missing";
        print_help();
        return EXIT_FAILURE;
    }

    if ((gaf == true && binary_path.empty()) || (gaf == true && pg_path.empty())) {
        cerr << "GAF file can be generated only with binary phenotype AND with the pg graph";
        print_help();
        return EXIT_FAILURE;
    }

    if (phenotype > 1) {
        cerr << "Only one kind of analysis/phenotype is allowed";
        print_help();
        return EXIT_FAILURE;
    }

    // scope declaration
    std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> snarls_chr;
    std::unique_ptr<bdsg::PackedGraph> pg;

    if (!snarl_path.empty()){
        snarls_chr = parse_snarl_path(snarl_path);

    } else if (!pg_path.empty() && !dist_path.empty()) {
        std::cout << "Start snarl analysis... " << std::endl;
        auto start_0 = std::chrono::high_resolution_clock::now();
        auto [stree, pg, root, pp_overlay] = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
        string output_file = output_dir + "/snarl_analyse.tsv";
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, true);
        auto end_0 = std::chrono::high_resolution_clock::now();
        std::cout << "Snarl analysis : " << std::chrono::duration<double>(end_0 - start_0).count() << " s" << std::endl;
    }

    auto [list_samples, ptr_vcf, hdr, rec] = parseHeader(vcf_path);    
    std::unordered_map<std::string, bool> binary;
    std::unordered_map<std::string, double> quantitative;

    if (!binary_path.empty()) {
        check_format_binary_phenotype(binary_path);
        binary = parse_binary_pheno(binary_path);
        check_match_samples(binary, list_samples);

        string output_binary = output_dir + "/binary_gwas.tsv";
        std::ofstream outf(output_binary, std::ios::binary);
        std::string headers = "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n";
        outf.write(headers.c_str(), headers.size());

        chromosome_chuck_binary(ptr_vcf, hdr, rec, list_samples, snarls_chr, binary, outf);
        if (gaf) {
            string output_gaf = output_dir + "/snarl.gaf";
            parse_input_file(output_binary, snarls_chr, *pg, output_gaf);
        }

    } else if (!quantitative_path.empty()) {
        check_format_quantitative_phenotype(quantitative_path);
        quantitative = parse_quantitative_pheno(quantitative_path);
        check_match_samples(quantitative, list_samples);

        string quantitive_output = output_dir + "/quantitative_gwas.tsv";
        std::ofstream outf(quantitive_output, std::ios::binary);
        std::string headers = "CHR\tPOS\tSNARL\tTYPE\tBETA\tSE\tP\tALLELE_NUM\n";
        outf.write(headers.c_str(), headers.size());

        chromosome_chuck_quantitative(ptr_vcf, hdr, rec, list_samples, snarls_chr, quantitative, outf);

    } else if (!eqtl_path.empty()) {
        // check_format_eqtl_phenotype(eqtl_path);
        auto eqtl = parseEQTLFile(eqtl_path);
        // check_match_samples(eqtl, list_samples);

        // string eqtl_output = output_dir + "/eqtl_gwas.tsv";
        // std::ofstream outf(eqtl_output, std::ios::binary);
        // std::string headers = "CHR\tPOS\tSNARL\tTYPE\tSE\tBETA\tP\n";
        // outf.write(headers.c_str(), headers.size());

        // chromosome_chuck_quantitative(ptr_vcf, hdr, rec, list_samples, snarls_chr, eqtl, outf);
    }

    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Time Gwas analysis : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;

    return EXIT_SUCCESS;
}

// BINARY
// ./stoat_cxx -p ../data/binary/pg.pg -d ../data/binary/pg.dist -v ../data/binary/binary.vcf.gz -b ../data/binary/phenotype.tsv

// QUANTITATIVE
// ./stoat_cxx -p ../data/quantitative/pg.pg -d ../data/quantitative/pg.dist -v ../data/quantitative/quantitative.vcf.gz -b ../data/quantitative/phenotype.tsv
