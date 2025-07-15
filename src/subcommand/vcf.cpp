// This file is part of STOAT 0.0.1, copyright (C) 2024-2025 
// Authors : Matis Alias-Bagarre, Jean Monlong & Xian-hui Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

#include "../snarl_data_t.hpp"
#include "../snarl_analyzer.hpp"
#include "../arg_parser.hpp"
#include "../matrix.hpp"
#include "../gaf_creator.hpp"
#include "../post_processing.hpp"

namespace stoat_vcf {

void print_help_vcf() {
    std::cerr << "Usage: stoat vcf [options]\n\n"
              << "  -p, --pg FILE                Path to the packed graph file (.pg)\n"
              << "  -d, --dist FILE              Path to the packed distance index file (.dist)\n"
              << "  -v, --vcf FILE               Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -s, --snarl FILE             Path to the snarl file (.txt or .tsv)\n"
              << "  -r, --chr FILE               Path to the chromosome reference file (.txt)\n"
              << "  -b, --binary FILE            Path to the binary phenotype group file (.txt or .tsv)\n"
              << "  -q, --quantitative FILE      Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  -e, --eqtl FILE              Path to the Expression Quantitative Trait Loci file (.txt or .tsv)\n"
              << "  -m, --make-bed                   Create plink format files (.bed, .bim, .fam)\n"
              << "  -c, --covariate FILE             Path to the covariate file (.txt or .tsv)\n"
              << "  -C, --covar-name NAME            Covariate column name(s) used for GWAS (comma-separated if multiple)\n"
              << "  -k, --kinship FILE           Path to the kinship matrix file (.txt or .tsv)\n"
              << "  -g, --gaf                    Generate a GAF file from GWAS results\n"
              << "  -i, --children INT               Max number of children per snarl in decomposition (default: 50)\n"
              << "  -y, --cycle INT                  Max number of authorized cycles in snarl decomposition (default: 1)\n"
              << "  -l, --path-length INT            Max number of nodes in paths during snarl decomposition (default: 10,000)\n"
              << "  -G, --gene-position FILE     Path to the gene position file (.txt or .tsv)\n"
              << "  -w, --windows-gene INT       Window length from gene boundaries for snarl inclusion in eQTL (default: 1,000,000)\n"
              << "  -T, --table-threshold FLOAT  P-value threshold for regression table output (default: disabled)\n"
              << "  -M, --maf FLOAT                  Minimum allele frequency threshold (default: 0.01)\n"
              << "  -t, --thread INT             Number of threads to use (default: 1)\n"
              << "  -o, --output DIR             Output directory name (VCF GWAS mode)\n"
              << "  -h, --help                   Print this help message\n";
}

int stoat_vcf(int argc, char* argv[]) {
    
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, pg_path, dist_path, 
        chromosome_path, binary_path, quantitative_path, 
        eqtl_path, covariate_path, gene_position_path, 
        kinship_path, output_dir;

    size_t phenotype = 0;
    size_t cycle_threshold = 1;
    size_t children_threshold = 50;
    size_t path_length_threshold = 10000;
    size_t windows_gene_threshold = 1000000;
    double table_threshold = -1;
    double maf_threshold = 0.05;
    bool gaf = false;
    bool only_snarl_parsing = false;
    bool show_help = false;
    bool make_bed = false;

    std::vector<std::string> covar_names;

    // Parse arguments
    int c;

    static struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"snarl", required_argument, 0, 's'},
        {"pg", required_argument, 0, 'p'},
        {"dist", required_argument, 0, 'd'},
        {"chr", required_argument, 0, 'r'},
        {"binary", required_argument, 0, 'b'},
        {"quantitative", required_argument, 0, 'q'},
        {"eqtl", required_argument, 0, 'e'},
        {"make-bed", no_argument, 0, 'm'},
        {"covariate", required_argument, 0, 'c'},
        {"covar-name", required_argument, 0, 'C'},
        {"kinship", required_argument, 0, 'k'},
        {"gaf", no_argument, 0, 'g'},
        {"children", required_argument, 0, 'i'},
        {"cycle", required_argument, 0, 'y'},
        {"path-length", required_argument, 0, 'l'},
        {"gene-position", required_argument, 0, 'G'},
        {"windows-gene", required_argument, 0, 'w'},
        {"table-threshold", required_argument, 0, 'T'},
        {"maf", required_argument, 0, 'M'},
        {"thread", required_argument, 0, 't'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "v:s:p:d:r:b:q:e:m:c:C:k:g:i:y:l:G:w:T:M:t:o:h", long_options, nullptr)) != -1) {
        switch (c) {
            case 'v': vcf_path = optarg; stoat_vcf::check_file(vcf_path); break;
            case 's': snarl_path = optarg; stoat_vcf::check_file(snarl_path); break;
            case 'p': pg_path = optarg; stoat_vcf::check_file(pg_path); break;
            case 'd': dist_path = optarg; stoat_vcf::check_file(dist_path); break;
            case 'r': chromosome_path = optarg; stoat_vcf::check_file(chromosome_path); break;
            case 'b': binary_path = optarg; phenotype++; stoat_vcf::check_file(binary_path); break;
            case 'q': quantitative_path = optarg; phenotype++; stoat_vcf::check_file(quantitative_path); break;
            case 'e': eqtl_path = optarg; phenotype++; stoat_vcf::check_file(eqtl_path); break;
            case 'm': make_bed = true; break;
            case 'c': covariate_path = optarg; stoat_vcf::check_file(covariate_path); break;
            case 'C': {
                std::stringstream ss(optarg);
                std::string token;
                while (std::getline(ss, token, ',')) covar_names.push_back(token);
                break;
            }
            case 'k': kinship_path = optarg; stoat_vcf::check_file(kinship_path); break;
            case 'g': gaf = true; break;
            case 'i':
                children_threshold = std::stoi(optarg);
                if (children_threshold < 2) {
                    std::cerr << "Error: Children threshold must be > 1\n";
                    return EXIT_FAILURE;
                }
                break;
            case 'y':
                cycle_threshold = std::stoi(optarg);
                if (cycle_threshold < 1) {
                    std::cerr << "Error: Cycle threshold must be > 0\n";
                    return EXIT_FAILURE;
                }
                break;
            case 'l':
                path_length_threshold = std::stoi(optarg);
                if (path_length_threshold < 2) {
                    std::cerr << "Error: Path length threshold must be > 1\n";
                    return EXIT_FAILURE;
                }
                break;
            case 'G': gene_position_path = optarg; stoat_vcf::check_file(gene_position_path); break;
            case 'w':
                windows_gene_threshold = std::stoi(optarg);
                if (windows_gene_threshold < 1) {
                    std::cerr << "Error: Windows gene threshold must be > 0\n";
                    return EXIT_FAILURE;
                }
                break;
            case 'T':
                table_threshold = std::stod(optarg);
                if (table_threshold <= 0 || table_threshold > 1) {
                    std::cerr << "Error: Table threshold must be in (0,1]\n";
                    return EXIT_FAILURE;
                }
                break;
            case 'M':
                maf_threshold = std::stod(optarg);
                if (maf_threshold < 0 || maf_threshold > 1) {
                    std::cerr << "Error: MAF must be in [0,1]\n";
                    return EXIT_FAILURE;
                }
                break;
            case 't':
                if (std::stoi(optarg) < 1) {
                    std::cerr << "Error: Number of threads must be > 0\n";
                    return EXIT_FAILURE;
                }
                omp_set_num_threads(std::stoi(optarg));
                break;
            case 'o': output_dir = optarg; break;
            case 'h': print_help_vcf(); exit(EXIT_SUCCESS); break;
            default:
                std::cerr << "Unknown argument. Use -h or --help for usage.\n";
                return EXIT_FAILURE;
        }
    }

    if (show_help || argc == 2) {
        print_help_vcf();
        return EXIT_FAILURE;
    }
    
    if (output_dir.empty()) {
        output_dir = "output";
    }

    if (!covariate_path.empty() && covar_names.empty()) {
        std::cerr << "If --covariate path is provided you must add the column name(s), using --covar-name" << "\n";
        print_help_vcf();
        return EXIT_FAILURE;
    }

    if ((!eqtl_path.empty() && gene_position_path.empty()) || (eqtl_path.empty() && !gene_position_path.empty())) {
        std::cerr << "eqtl phenotype file and gene position file must be provided together" << "\n";
        print_help_vcf();
        return EXIT_FAILURE;
    }

    auto start_1 = std::chrono::high_resolution_clock::now();
    std::filesystem::create_directory(output_dir);
    
    if (chromosome_path.empty() && snarl_path.empty()) {
        std::cout << "Warning : chromosome_path file not provided, 'ref' reference chromosome name will be used instead" << std::endl;
    }

    std::unordered_set<std::string> ref_chr = (!chromosome_path.empty()) ? stoat_vcf::parse_chromosome_reference(chromosome_path) : std::unordered_set<std::string>{"ref"};
    std::string regression_dir = output_dir + "/regression";

    if (table_threshold != -1) {
        std::filesystem::create_directory(regression_dir);
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
        std::cerr << "There are 5 ways to lauch stoat : " << std::endl;
        std::cerr << "Case 1: snarl_path + vcf_path + phenotype (+ optional file)" << std::endl;
        std::cerr << "Case 2: pg_path + dist_path + vcf_path + phenotype (+ optional file)" << std::endl;
        std::cerr << "Case 3: pg_path + dist_path" << std::endl;
        std::cerr << "Case 4: pg_path + dist_path + vcf_path + --make-bed" << std::endl;
        std::cerr << "Case 5: snarl_path + vcf_path + --make-bed" << std::endl;

        print_help_vcf();
        return EXIT_FAILURE;
    }

    if ((gaf == true && binary_path.empty()) || (gaf == true && pg_path.empty())) {
        std::cerr << "GAF file can be generated only with binary phenotype AND with the pg graph" << std::endl;
        print_help_vcf();
        return EXIT_FAILURE;
    }

    std::vector<std::string> list_samples;
    htsFile* ptr_vcf;
    bcf_hdr_t* hdr;
    bcf1_t* rec;

    if (!only_snarl_parsing) {
        std::tie(list_samples, ptr_vcf, hdr, rec) = stoat_vcf::parseHeader(vcf_path); 
    }

    //////////////////// Load the phenotypes and covariate matrix from files

    std::vector<bool> binary_phenotype;
    std::vector<double> quantitative_phenotype;

    // dict chr:string : vector{geneName:string, sample_expression:vector<double>, start_pos:size_t, end_pos:size_t}
    std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>> eqtl_phenotype;
    std::vector<std::vector<double>> covariate;

    if (!covariate_path.empty()) {
        covariate = stoat_vcf::parse_covariates(covariate_path, covar_names, list_samples);
        covariate = stoat_vcf::parse_covariates(covariate_path, covar_names, list_samples);
    }

    if (!binary_path.empty()) {
        binary_phenotype = stoat_vcf::parse_binary_pheno(binary_path, list_samples);

    } else if (!quantitative_path.empty()) {
        quantitative_phenotype = stoat_vcf::parse_quantitative_pheno(quantitative_path, list_samples);

    } else if (!eqtl_path.empty() && !gene_position_path.empty()) {
        eqtl_phenotype = stoat_vcf::parse_qtl_gene_file(eqtl_path, gene_position_path, list_samples);
    }

    stoat_vcf::KinshipMatrix kinship;
    if (!kinship_path.empty()) {
        kinship.parseKinshipMatrix(kinship_path);
    }

    // Load or calculate the snarl information
    // scope declaration
    // chr : <snarl, paths, pos(start, end), type>
    std::unordered_map<std::string, std::vector<stoat_vcf::Snarl_data_t>> snarls_chr;
    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<bdsg::PackedGraph> pg;
    handlegraph::net_handle_t root;
    std::unique_ptr<bdsg::PackedPositionOverlay> pp_overlay;

    if (!snarl_path.empty()){ // If we have already saved the paths in snarls, load them
        snarls_chr = stoat_vcf::parse_snarl_path(snarl_path);
    } else { // Otherwise, find them from the graph and snarl tree
        std::cout << "Start snarl analysis... " << std::endl;
        auto start_0 = std::chrono::high_resolution_clock::now();
        // Load the snarl tree and graph
        std::tie(stree, pg, root, pp_overlay) = stoat_vcf::parse_graph_tree(pg_path, dist_path);

        // std::vector<std::tuple<handlegraph::net_handle_t, std::string, size_t, size_t, bool>>
        // snarl_net_grah, chr_ref, start_pos, end_pos, is_on_ref
        auto snarls = stoat_vcf::save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);

        std::string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
        std::string output_file = output_dir + "/snarl_analyse.tsv";

        // Go through snarls and fill in snarls_chr 
        snarls_chr = stoat_vcf::loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, cycle_threshold, only_snarl_parsing);
        auto end_0 = std::chrono::high_resolution_clock::now();
        std::cout << "Snarl decomposition : " << std::chrono::duration<double>(end_0 - start_0).count() << " s" << std::endl;
        if (only_snarl_parsing) {
            return EXIT_SUCCESS;
        }

        // Clean up unique_ptr except pg
        stree.reset();
        pp_overlay.reset();
    }

    //////////////////////////////////////// Go through the vcf, do the analysis, and write the output

    auto start_2 = std::chrono::high_resolution_clock::now();

    std::shared_ptr<stoat_vcf::SnarlAnalyzer> snarl_analyzer;
    stoat_vcf::EdgeBySampleMatrix edge_matrix_empty(list_samples, 0, 0);
    stoat::phenotype_type_t phenotype_type;

    // Decide which type of SnarlAnalyzer we want
    if (!binary_path.empty()) {
        // binary
        if (!covariate.empty()){
            // Binary covariate
            snarl_analyzer.reset(new stoat_vcf::BinaryCovarSnarlAnalyzer(snarls_chr, list_samples, covariate, maf_threshold, table_threshold, binary_phenotype));
        } else {
            // Binary normal
            snarl_analyzer.reset(new stoat_vcf::BinarySnarlAnalyzer(snarls_chr, list_samples, maf_threshold, table_threshold, binary_phenotype));
        }
        phenotype_type = stoat::BINARY; 
    } else if (!quantitative_path.empty()) {
        // Quantitative
        snarl_analyzer.reset(new stoat_vcf::QuantitativeSnarlAnalyzer(snarls_chr, list_samples, covariate, maf_threshold, table_threshold, quantitative_phenotype));
        phenotype_type = stoat::QUANTITATIVE; 
    } else if (!eqtl_path.empty()) {
        // EQTL
        snarl_analyzer.reset(new stoat_vcf::EQTLSnarlAnalyzer(snarls_chr, list_samples, covariate, maf_threshold, table_threshold, eqtl_phenotype, windows_gene_threshold));
        phenotype_type = stoat::EQTL; 
    }

    std::string output_tsv = output_dir + (phenotype_type == stoat::BINARY       ? "/binary_table.tsv" : 
                                            (phenotype_type == stoat::QUANTITATIVE ? "/quantitative_table.tsv" 
                                                                                        : "/eqtl_gwas.tsv"));

    snarl_analyzer->process_snarls_by_chromosome_chunk(ptr_vcf, hdr, rec, edge_matrix_empty, regression_dir, output_tsv);

    std::string output_significative = output_dir + (phenotype_type == stoat::BINARY       ?  "/top_variant_binary.tsv" : 
                                                    (phenotype_type == stoat::QUANTITATIVE ? "/top_variant_quantitative.tsv" 
                                                                                                : "/top_variant_eqtl.tsv"));

    stoat_vcf::add_BH_adjusted_column(output_tsv, output_dir, output_significative, phenotype_type);

    if (phenotype_type == stoat::BINARY && gaf) {
        std::string output_gaf = output_dir + "/binary_table.gaf";
        stoat_vcf::gaf_creation(output_tsv, snarls_chr, *pg, output_gaf);
    }

    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Snarl analysis : " << std::chrono::duration<double>(end_1 - start_2).count() << " s" << std::endl;
    std::cout << "Time Gwas analysis : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;
    return EXIT_SUCCESS;
}

} // end stoat_vcf
