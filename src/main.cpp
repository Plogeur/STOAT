// This file is part of STOAT 0.0.1, copyright (C) 2024-2025 Matis Alias-Bagarre, Jean Monlong & Xian-hui Chang.
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
#include <chrono>
#include <Eigen/Dense>
#include <cstdlib>
#include <getopt.h>

#include "arg_parser.hpp"
#include "snarl_parser.hpp"     
#include "matrix.hpp"
#include "list_snarl_paths.hpp"
#include "gaf_creator.hpp"
#include "post_processing.hpp"

using namespace std;

// 
// stoat vcf ...
// stoat graph ... 
void print_help() {
    std::cout << "Usage: stoat [options]\n\n"
              << "  -p, --pg FILE                Path to the packed graph file (.pg)\n"
              << "  -d, --dist FILE              Path to the packed distance index file (.dist)\n"
              << "  -v, --vcf FILE               Path to the VCF file (.vcf or .vcf.gz)\n"
              << "  -s, --snarl FILE             Path to the snarl file (.txt or .tsv)\n"
              << "  -r, --chr FILE               Path to the chromosome reference file (.txt)\n"
              << "  -b, --binary FILE            Path to the binary phenotype group file (.txt or .tsv)\n"
              << "  -q, --quantitative FILE      Path to the quantitative phenotype file (.txt or .tsv)\n"
              << "  -e, --eqtl FILE              Path to the Expression Quantitative Trait Loci file (.txt or .tsv)\n"
              << "  -m, --make-bed               Create plink format files (.bed, .bim, .fam)\n"
              << "  -c, --covariate FILE         Path to the covariate file (.txt or .tsv)\n"
              << "  -C, --covar-name NAME        Covariate column name(s) used for GWAS (comma-separated if multiple)\n"
              << "  -k, --kinship FILE           Path to the kinship matrix file (.txt or .tsv)\n"
              << "  -g, --gaf                    Generate a GAF file from GWAS results\n"
              << "  -i, --children INT           Max number of children per snarl in decomposition (default: 50)\n"
              << "  -y, --cycle INT              Max number of authorized cycles in snarl decomposition (default: 1)\n"
              << "  -l, --path-length INT        Max number of nodes in paths during snarl decomposition (default: 10,000)\n"
              << "  -G, --gene-position FILE     Path to the gene position file (.txt or .tsv)\n"
              << "  -w, --windows-gene INT       Window length from gene boundaries for snarl inclusion in eQTL (default: 1,000,000)\n"
              << "  -T, --table-threshold FLOAT  P-value threshold for regression table output (default: disabled)\n"
              << "  -M, --maf FLOAT              Minimum allele frequency threshold (default: 0.01)\n"
              << "  -t, --thread INT             Number of threads to use (default: 1)\n"
              << "  -o, --output DIR             Output directory name (VCF GWAS mode)\n"
              << "  -h, --help                   Print this help message\n";
}


int main(int argc, char* argv[]) {
    // Declare variables to hold argument values
    std::string vcf_path, snarl_path, pg_path, dist_path, 
        chromosome_path, binary_path, quantitative_path, 
        eqtl_path, covariate_path, gene_position_path, 
        kinship_path, output_dir;

    size_t num_threads = 1;
    size_t phenotype = 0;
    size_t cycle_threshold = 1;
    size_t children_threshold = 50;
    size_t path_length_threshold = 10000;
    size_t windows_gene_threshold = 1000000;
    double table_threshold = -1;
    double maf = 0.99; // inversed MAF 
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
            case 'v': vcf_path = optarg; check_file(vcf_path); break;
            case 's': snarl_path = optarg; check_file(snarl_path); break;
            case 'p': pg_path = optarg; check_file(pg_path); break;
            case 'd': dist_path = optarg; check_file(dist_path); break;
            case 'r': chromosome_path = optarg; check_file(chromosome_path); break;
            case 'b': binary_path = optarg; phenotype++; check_file(binary_path); break;
            case 'q': quantitative_path = optarg; phenotype++; check_file(quantitative_path); break;
            case 'e': eqtl_path = optarg; phenotype++; check_file(eqtl_path); break;
            case 'm': make_bed = true; break;
            case 'c': covariate_path = optarg; check_file(covariate_path); break;
            case 'C': {
                std::stringstream ss(optarg);
                std::string token;
                while (std::getline(ss, token, ',')) covar_names.push_back(token);
                break;
            }
            case 'k': kinship_path = optarg; check_file(kinship_path); break;
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
            case 'G': gene_position_path = optarg; check_file(gene_position_path); break;
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
                maf = 1 - std::stod(optarg);
                if (maf < 0 || maf > 1) {
                    std::cerr << "Error: MAF must be in [0,1]\n";
                    return EXIT_FAILURE;
                }
                break;
            case 't':
                num_threads = std::stoi(optarg);
                if (num_threads < 1) {
                    std::cerr << "Error: Number of threads must be > 0\n";
                    return EXIT_FAILURE;
                }
                break;
            case 'o': output_dir = optarg; break;
            case 'h': print_help(); exit(EXIT_SUCCESS); break;
            default:
                std::cerr << "Unknown argument. Use -h or --help for usage.\n";
                return EXIT_FAILURE;
        }
    }

    if (show_help || argc == 1) {
        print_help();
        return EXIT_FAILURE;
    }
    
    if (output_dir.empty()) {
        output_dir = "output";
    }

    if (!covariate_path.empty() && covar_names.empty()) {
        std::cerr << "If --covariate path is provided you must add the column name(s), using --covar-name" << "\n";
        print_help();
        return EXIT_FAILURE;
    }

    if ((!eqtl_path.empty() && gene_position_path.empty()) || (eqtl_path.empty() && !gene_position_path.empty())) {
        std::cerr << "eqtl phenotype file and gene position file must be provided together" << "\n";
        print_help();
        return EXIT_FAILURE;
    }

    auto start_1 = std::chrono::high_resolution_clock::now();
    std::filesystem::create_directory(output_dir);
    
    if (chromosome_path.empty() && snarl_path.empty()) {
        std::cout << "Warning : chromosome_path file not provided, 'ref' reference chromosome name will be used instead" << std::endl;
    }

    std::unordered_set<std::string> ref_chr = (!chromosome_path.empty()) ? parse_chromosome_reference(chromosome_path) : std::unordered_set<std::string>{"ref"};
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
        std::cerr << "GAF file can be generated only with binary phenotype AND with the pg graph" << std::endl;
        print_help();
        return EXIT_FAILURE;
    }

    std::vector<std::string> list_samples;
    htsFile* ptr_vcf;
    bcf_hdr_t* hdr;
    bcf1_t* rec;

    if (!only_snarl_parsing) {
        std::tie(list_samples, ptr_vcf, hdr, rec) = parseHeader(vcf_path); 
    }

    std::vector<bool> binary;
    std::vector<double> quantitative;
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>> eqtl;
    std::vector<std::vector<double>> covariate;

    if (!covariate_path.empty()) {
        covariate = parse_covariates(covariate_path, covar_names, list_samples);
    }

    if (!binary_path.empty()) {
        binary = parse_binary_pheno(binary_path, list_samples);

    } else if (!quantitative_path.empty()) {
        quantitative = parse_quantitative_pheno(quantitative_path, list_samples);

    } else if (!eqtl_path.empty() && !gene_position_path.empty()) {
        eqtl = parse_qtl_gene_file(eqtl_path, gene_position_path, list_samples);
    }

    KinshipMatrix kinship;
    if (!kinship_path.empty()) {
        // check_format_kinship(kinship_path);
        kinship = parseKinshipMatrix(kinship_path);
    }

    // scope declaration
    // chr : <snarl, paths, pos(start, end), type>
    // TODO : replace std::tuple<string, vector<string>, size_t, size_t, vector<string>> to 5 vector (to reduce space/time usage)
    std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>> snarls_chr;
    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<bdsg::PackedGraph> pg;
    handlegraph::net_handle_t root;
    std::unique_ptr<bdsg::PackedPositionOverlay> pp_overlay;

    if (!snarl_path.empty()){
        snarls_chr = parse_snarl_path(snarl_path);
    } else {

        std::cout << "Start snarl analysis... " << std::endl;
        auto start_0 = std::chrono::high_resolution_clock::now();
        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);

        // vector<tuple<net_handle_t, string, size_t, size_t, bool>>
        // snarl_net_grah, chr_ref, start_pos, end_pos, is_on_ref
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);

        string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
        string output_file = output_dir + "/snarl_analyse.tsv";

        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, cycle_threshold, only_snarl_parsing);
        auto end_0 = std::chrono::high_resolution_clock::now();
        std::cout << "Snarl decomposition : " << std::chrono::duration<double>(end_0 - start_0).count() << " s" << std::endl;
        if (only_snarl_parsing) {
            return EXIT_SUCCESS;
        }

        // Clean up unique_ptr except pg
        stree.reset();
        pp_overlay.reset();
    }

    auto start_2 = std::chrono::high_resolution_clock::now();

    if (make_bed) {

        std::vector<std::pair<std::string, int>> pheno;
        for (const auto& sample : list_samples) {
            pheno.push_back({sample, -9}); // initilize all phenotypes to -9
        }

        const std::string output_fam = output_dir + ".fam";
        create_fam(pheno, output_fam);
        chromosome_chuck_make_bed(ptr_vcf, hdr, rec, list_samples, snarls_chr, output_dir);

        auto end_1 = std::chrono::high_resolution_clock::now();
        std::cout << "Time genotype plink files creations : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;
        return EXIT_SUCCESS;

    } else if (!binary_path.empty()) {

        string output_binary = output_dir + "/binary_analysis.tsv";
        chromosome_chuck_binary(ptr_vcf, hdr, rec, list_samples, snarls_chr, binary, covariate, maf, kinship, num_threads, table_threshold, regression_dir, output_binary);

        string output_significative = output_dir + "/top_variant_binary.tsv";
        string phenotype_type = covariate.empty() ? "binary" : "quantitative";
        add_BH_adjusted_column(output_binary, output_significative, phenotype_type);

        if (gaf) {
            string output_gaf = output_dir + "/binary_analysis.gaf";
            gaf_creation(output_binary, snarls_chr, *pg, output_gaf);
        }

    } else if (!quantitative_path.empty()) {

        string output_quantitive = output_dir + "/quantitative_analysis.tsv";
        chromosome_chuck_quantitative(ptr_vcf, hdr, rec, list_samples, snarls_chr, quantitative, covariate, maf, kinship, num_threads, table_threshold, regression_dir, output_quantitive);

        string output_significative = output_dir + "/top_variant_quantitative.tsv";
        string phenotype_type = "quantitative";
        add_BH_adjusted_column(output_quantitive, output_significative, phenotype_type);

    } else if (!eqtl_path.empty()) {

        string eqtl_output = output_dir + "/eqtl_gwas.tsv";
        chromosome_chuck_eqtl(ptr_vcf, hdr, rec, list_samples, snarls_chr, eqtl, covariate, maf, 
            kinship, num_threads, table_threshold, regression_dir, windows_gene_threshold, eqtl_output);
        
        string output_significative = output_dir + "/top_variant_eqtl.tsv";
        string phenotype_type = "eqtl";
        add_BH_adjusted_column(eqtl_output, output_significative, phenotype_type);
    }

    auto end_1 = std::chrono::high_resolution_clock::now();
    std::cout << "Snarl analysis : " << std::chrono::duration<double>(end_1 - start_2).count() << " s" << std::endl;
    std::cout << "Time Gwas analysis : " << std::chrono::duration<double>(end_1 - start_1).count() << " s" << std::endl;
    return EXIT_SUCCESS;
}

// DROSO
// ./stoat_cxx -p ../data_droso/fly.pg -d ../data_droso/fly.dist -v ../data_droso/merged.vcf -q ../data_droso/phenotype.tsv --output ../output_droso
   
// DROSO
// ./stoat_cxx -p ../data/droso/fly.pg -d ../data/droso/fly.dist -r ../data/droso/chromosome_ref.tsv --output ../output_droso
// sed -i 's/dm6#0#chr2L/1/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr2R/2/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr3L/3/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr3R/4/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chr4/5/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrX/6/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrY/7/g' ../output_droso/snarl_analyse.tsv
// sed -i 's/dm6#0#chrM/8/g' ../output_droso/snarl_analyse.tsv
// ./stoat_cxx -s ../output_droso/snarl_analyse.tsv -v ../data/droso/merging_stoat.vcf -q ../data/droso/pangenome_pheno.tsv --output ../output_droso

// BINARY
// ./stoat_cxx -p ../data/binary/pg.pg -d ../data/binary/pg.dist -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output

// BINARY + COVARIATE
// ./stoat_cxx -p ../data/binary/pg.pg -d ../data/binary/pg.dist -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --covariate ../data/binary/covariate.tsv --covar-name CP1,SEX,CP3 --output ../output

// QUANTITATIVE
// ./stoat_cxx -p ../data/quantitative/pg.pg -d ../data/quantitative/pg.dist -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv --output ../output

// QUANTITATIVE + COVARIATE
// ./stoat_cxx -p ../data/quantitative/pg.pg -d ../data/quantitative/pg.dist -v ../data/quantitative/merged_output.vcf.gz -q ../data/quantitative/phenotype.tsv  --covariate ../data/quantitative/covariate.tsv --covar-name CP1,SEX,CP3 --output ../output

// EQTL
// ./stoat_cxx -s ../test_data/quantitative/paths_snarl.tsv -v ../test_data/quantitative/variants.vcf -e ../test_data/quantitative/qtl.tsv --gene-position ../test_data/quantitative/gene_position.tsv --output ../output

// TEST
// ./stoat_cxx -p ../tests/graph_test/3th_snp.pg -d ../tests/graph_test/3th_snp.dist --output ../output

// BINARY-PLINK
// ./stoat_cxx -p ../data/binary/pg.pg -d ../data/binary/pg.dist -v ../data/binary/merged_output.vcf.gz --make-bed --output ../output

// QUANTITATIVE-PLINK
// ./stoat_cxx -p ../data/quantitative/pg.pg -d ../data/quantitative/pg.dist -v ../data/quantitative/merged_output.vcf.gz --make-bed --output ../output

// SIMULATION NEW
// ./stoat_cxx -v ../data/simu/variants.vcf -s ../data/simu/paths_snarl.tsv -b ../data/simu/phenotypes.txt --covariate ../data/simu/covar.tsv --covar-name AGE,SEX,PC1,PC2 --output ../output

// ./stoat_cxx -v ../data/simu/variants.vcf -s ../data/simu/paths_snarl.tsv -b ../data/simu/phenotypes.txt --make-bed --output ../output
// plink --bfile ../output/output --pheno ../data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/stoat_plink

// PLINK
// plink --vcf ../data/simu/variants.vcf --make-bed --allow-extra-chr --out ../output/genotype
// plink --bfile ../output/genotype --pheno ../data/simu/phenotypes.txt --pheno-name PHENO --assoc --allow-no-sex --allow-extra-chr --out ../output/plink

// VALGRIND
// valgrind --tool=callgrind ./stoat_cxx -s ../data/binary/snarl_paths.tsv -v ../data/binary/merged_output.vcf.gz -b ../data/binary/phenotype.tsv --output ../output
// kcachegrind callgrind.out.<id>