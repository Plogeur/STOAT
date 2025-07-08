#include "snarl_analyzer.hpp"
#include "matrix.hpp"
#include "binary_table.hpp"
#include "quantitative_table.hpp"
#include "utils.hpp"
#include "arg_parser.hpp"
#include "writer.hpp"
#include "omp.h"

namespace stoat_vcf {

SnarlAnalyzer::SnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, double maf, double table_threshold) :
chr_to_snarl_data(chr_to_snarl_data), list_samples(list_samples), covariate(covariate), maf(maf), table_threshold(table_threshold), edge_matrix(list_samples,0,0) {};

BinarySnarlAnalyzer::BinarySnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
    double maf, double table_threshold, const std::vector<bool>& binary_phenotype) :
SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf, table_threshold), binary_phenotype(binary_phenotype) {};

BinaryCovarSnarlAnalyzer::BinaryCovarSnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, const std::vector<bool>& binary_phenotype) :
SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf, table_threshold), binary_phenotype(binary_phenotype) {};

QuantitativeCovarSnarlAnalyzer::QuantitativeSnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, const std::vector<double>& quantitative_phenotype) :
SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf, table_threshold), quantitative_phenotype(quantitative_phenotype) {};

EQTLCovarSnarlAnalyzer::EQTLSnarlAnalyzer(const std::unordered_map<std::string, std::vector<Snarl_data_t>>& chr_to_snarl_data, const std::vector<std::string>& list_samples, 
    const std::vector<std::vector<double>>& covariate, double maf, double table_threshold, 
    const std::unordered_map<std::string, std::vector<stoat_vcf::Qtl_data>>& eqtl_map,
    size_t windows_gene_threshold) :
SnarlAnalyzer(chr_to_snarl_data, list_samples, covariate, maf, table_threshold), eqtl_map(eqtl_map) {};

void BinarySnarlAnalyzer::write_header(std::ofstream& outf) {
    write_binary_header(outf);
}
void BinaryCovarSnarlAnalyzer::write_header(std::ofstream& outf) {
    write_binary_covar_header(outf);
}
void QuantitativeSnarlAnalyzer::write_header(std::ofstream& outf) {
    write_quantitative_header(outf);
}
void EQTLSnarlAnalyzer::write_header(std::ofstream& outf) {
    write_eqtl_header(outf);
}

void SnarlAnalyzer::process_snarls_by_chromosome_chunk(
    htsFile* &ptr_vcf,
    bcf_hdr_t* &hdr,
    bcf1_t* &rec,
    const std::string& regression_dir_,
    const std::string& output_filename) {
    
    regression_dir = regression_dir_;
    outf(output_filename, std::ios::binary);

    // Write the header
    write_header(outf);

    // Go through the vcf and get chunks by chromosome. 
    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        chr = bcf_hdr_id2name(hdr, rec->rid);
        // Skip chromosomes not in chr_to_snarl_data
        while (chr_to_snarl_data.find(chr) == chr_to_snarl_data.end()) {
            std::cerr << "Warning: Chromosome " << chr << " not found in snarl paths file. Skipping." << std::endl;

            bool found_new_chr = false;
            while (bcf_read(ptr_vcf, hdr, rec) >= 0) {
                std::string chr_next = bcf_hdr_id2name(hdr, rec->rid);
                if (chr_next != chr) {
                    chr = chr_next;  // Update to the new chromosome
                    found_new_chr = true;
                    break;
                }
            }

            if (!found_new_chr) {
                return;  // exit if no more records are available
            }
        }

        std::cout << "> " << chr << std::endl;
        size_t size_chr = chr_to_snarl_data.at(chr).size();

        // Make genotype matrix by chromosome    
        auto [ptr_vcf_new, hdr_new, rec_new] = make_edge_matrix(ptr_vcf, hdr, rec, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        const auto& snarls = chr_to_snarl_data.at(chr);

        #pragma omp parallel for schedule(static)

        // Make the snarl test analysis
        // Iterate over each snarl
        for (const Snarl_data_t& snarl_data_s : snarls) {
            analyze_and_write_snarl(snarl_data_s);
        }
    }

    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> SnarlAnalyzer::make_edge_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, std::string &chr, size_t &num_paths_chr) {

    edge_matrix.reset(list_samples, num_paths_chr*4, list_samples.size() * 2);

    // loop over the VCF file for each line and stop where chr is different
    do {
        bcf_unpack(rec, BCF_UN_STR);

        // Check the INFO field for LV (Level Variant) and skip if LV != 0
        int32_t *lv = nullptr;
        int n_lv = 0;

        // Extract LV field from INFO skip if variant is lv != 0 to avoid duplication paths/snarl variant analysis
        if (bcf_get_info_int32(hdr, rec, "LV", &lv, &n_lv) > 0) {
            if (lv[0] != 0) {
                free(lv);
                continue;
            }
        }
        free(lv);

        // Extract genotypes (GT)
        int ngt = 0;
        int32_t *gt = nullptr;
        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt);

        // Extract AT field from INFO
        char *at_str = nullptr;
        int nat = 0;
        nat = bcf_get_info_string(hdr, rec, "AT", &at_str, &nat);
        std::vector<std::string> path_list;

        std::string at_value(at_str);  // Convert C-string to C++ std::string
        free(at_str);  // Free HTSlib-allocated memory

        // Split by comma
        std::stringstream ss(at_value);
        std::string item;
        while (std::getline(ss, item, ',')) {
            path_list.push_back(item);
        }

        // Decompose snarl paths [vector std::string] into [vector vector Edge_t]
        // paths : >123>213<234,>123<234,>123<234<345
        // list_paths_edge : [[Edge_t(123, 213), Edge_t(213, 234)], [...]]
        const std::vector<std::vector<stoat_vcf::Edge_t>> list_paths_edge = decompose_path_list_str(path_list);

        for (int i = 0; i < rec->n_sample; ++i) {
            int idex_path_allele_1 = bcf_gt_allele(gt[i * 2]);
            int idex_path_allele_2 = bcf_gt_allele(gt[i * 2 + 1]);
            size_t col_idx = i * 2;

            if (idex_path_allele_1 != -1) { // Handle missing genotypes
                for (const auto &edge_path_1 : list_paths_edge[idex_path_allele_1]) {
                    edge_matrix.push_matrix(edge_path_1, col_idx);
                }
            }

            if (idex_path_allele_2 != -1) { // Handle missing genotypes
                for (const auto &edge_path_2 : list_paths_edge[idex_path_allele_2]) {
                    edge_matrix.push_matrix(edge_path_2, col_idx + 1);
                }
            }
        }
        free(gt);

    } while ((bcf_read(ptr_vcf, hdr, rec) >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

    edge_matrix.shrink();
    return std::make_tuple(ptr_vcf, hdr, rec);
}

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>>& snarl_chr,
    const std::string& output_dir) {

    const std::string output_bed = output_dir + ".bed";
    const std::string output_bim = output_dir + ".bim";

    std::ofstream outbim(output_bim);
    std::ofstream outbed(output_bed, std::ios::binary);  // Open BED file as binary
    
    // Write the 3-byte 'BED' header for the BED file
    char bed_magic[] = {0x6C, 0x1B, 0x01};  // PLINK header: 0x6C ('l'), 0x1B, 0x01 (snp-major mode)
    outbed.write(bed_magic, 3);
    
    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        std::string chr = bcf_hdr_id2name(hdr, rec->rid);
        // Skip chromosomes not in snarl_chr
        while (snarl_chr.find(chr) == snarl_chr.end()) {
            std::cerr << "Warning: Chromosome " << chr << " not found in snarl paths file. Skipping." << std::endl;

            bool found_new_chr = false;
            while (bcf_read(ptr_vcf, hdr, rec) >= 0) {
                std::string chr_next = bcf_hdr_id2name(hdr, rec->rid);
                if (chr_next != chr) {
                    chr = chr_next;  // Update to the new chromosome
                    found_new_chr = true;
                    break;
                }
            }

            if (!found_new_chr) {
                return;  // exit if no more records are available
            }
        }

        std::cout << "> " << chr << std::endl;
        size_t size_chr = snarl_chr.at(chr).size();

        // Make genotype matrix by chromosome    
        auto [edge_matrix, ptr_vcf_new, hdr_new, rec_new] = make_edge_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr.at(chr);

        // Gwas analysis by chromosome
        create_bim_bed(snarl, list_samples.size(), edge_matrix, chr, outbim, outbed);
    }
    
    // Cleanup
    outbim.close();
    outbed.close();
    
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

std::pair<std::vector<size_t>, std::vector<size_t>> create_table_short_path(const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl, size_t sample_count,
    const EdgeBySampleMatrix& edge_matrix) {
    // TODO: Could just make the transposed version to start with?

    size_t length_column = list_path_snarl.size();
    std::vector<size_t> allele_number_list(length_column, 0);

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<size_t>> genotypes(sample_count, std::vector<size_t>(length_column, 0));

    // Genotype paths
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const stoat_vcf::Path_traversal_t& path_snarl = list_path_snarl[col_idx];
        std::vector<stoat_vcf::Edge_t> decomposed_snarl = decompose_path_to_edges(path_snarl);

        // Identify correct paths
        std::vector<size_t> idx_srr_save = identify_path(decomposed_snarl, edge_matrix, sample_count*2);

        for (size_t idx : idx_srr_save) {
            size_t srr_idx = idx / 2;  // Adjust index to correspond to the sample index
            genotypes[srr_idx][col_idx] += 1;
            allele_number_list[col_idx]++;
        }
    }

    std::vector<std::vector<size_t>> genotypes_transposed = transpose_matrix(genotypes);
    size_t major_index_1 = 0;
    size_t major_index_2 = 1;

    if (length_column > 2) {
        find_two_largest_indices(allele_number_list, major_index_1, major_index_2);
    }

    return {genotypes_transposed[major_index_1], genotypes_transposed[major_index_2]};
}

std::vector<std::vector<size_t>> transpose_matrix(const std::vector<std::vector<size_t>>& matrix) {
    if (matrix.empty()) return {};

    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    std::vector<std::vector<size_t>> transposed(cols, std::vector<size_t>(rows));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

void find_two_largest_indices(const std::vector<size_t>& vec, size_t& major_index_1, size_t& major_index_2) {

    // Ensure major_index_1 is the index of the larger of the first two
    if (vec[major_index_2] > vec[major_index_1]) {
        std::swap(major_index_1, major_index_2);
    }

    for (size_t i = 2; i < vec.size(); ++i) {
        if (vec[i] > vec[major_index_1]) {
            major_index_2 = major_index_1;
            major_index_1 = i;
        } else if (vec[i] > vec[major_index_2]) {
            major_index_2 = i;
        }
    }
}

void create_bim_bed(const std::vector<Snarl_data_t>& snarls, size_t sample_count, 
                    const EdgeBySampleMatrix& edge_matrix,
                    std::string chromosome, std::ofstream& outbim, std::ofstream& outbed) {

    // Iterate over each snarl
    // <snarl, paths, pos, type>
    for (const Snarl_data_t& snarl_data_s : snarls) {

        const std::vector<stoat_vcf::Path_traversal_t>& list_path_snarl = snarl_data_s.snarl_paths;
        size_t start_pos = snarl_data_s.start_positions;

        // if (list_path_snarl.size() > 2) {continue;} // avoid multiallelic var

        // Generate a genotype table for this snarl
        auto [allele_vector_0, allele_vector_1] = create_table_short_path(list_path_snarl, sample_count, edge_matrix);
        
        std::string allele1 = "A";  // Placeholder for allele 1
        std::string allele2 = "T";  // Placeholder for allele 2

        // chr id genetic_distance pos allele1 allele2
        outbim << chromosome << "\t" << pairToString(snarl_data_s.snarl_ids) << "\t0\t" << start_pos
                << "\t" << allele1 << "\t" << allele2 << "\n";
        
        // Write the genotypes for this SNP to the BED file
        unsigned char packed_byte = 0;  // A byte to store genotypes of 4 individuals
        int bit_pos = 0;

        // Loop through each sample (pair of alleles per individual)
        for (size_t snarl_list_idx = 0; snarl_list_idx < sample_count; ++snarl_list_idx) {

            size_t allele_0 = allele_vector_0[snarl_list_idx];      // number of allele for the individual for the first paths
            size_t allele_1 = allele_vector_1[snarl_list_idx];      // number of allele for the individual for the second paths

            // Encode the genotype as a 2-bit value based on the alleles
            unsigned char encoded_genotype = 0b10; // initialise as missing (./0, ./1, 1/. or 0/.)
            
            if (allele_0 == 1 || allele_1 == 1) {
                encoded_genotype = 0b01; // Heterozygous
            } else if (allele_0 == 2) {
                encoded_genotype = 0b00; // Homozygous major
            } else if (allele_1 == 2) {
                encoded_genotype = 0b11; // Homozygous minor
            }

            // Shift the encoded genotype into the correct position in the byte
            packed_byte |= (encoded_genotype << (bit_pos * 2));
            bit_pos++;

            // After 4 individuals, write the byte and reset the byte and bit position
            if (bit_pos == 4) {
                outbed.write(reinterpret_cast<char*>(&packed_byte), sizeof(unsigned char));
                packed_byte = 0;  // Reset the byte
                bit_pos = 0;      // Reset the bit position
            }
        }

        // If there are fewer than 4 individuals, write the remaining packed byte
        if (bit_pos > 0) {
            outbed.write(reinterpret_cast<char*>(&packed_byte), sizeof(unsigned char));
        }
    }
}

void create_fam(const std::vector<std::pair<std::string, int>> &pheno, 
    const std::string& output_path) {

    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        throw std::runtime_error("Unable to open output file: " + output_path);
    }

    for (const auto& [sample, phenotype] : pheno) {

        outfile << sample << " "           // FID (default to sample ID)
                << sample << " "           // IID (sample ID)
                << "0 0 0 "                // PID and MID (unknown)
                << phenotype << "\n";      // PHENOTYPE (-9 = missing)
    }
    outfile.close();
}

// Optimized function to extract node ID and update index
inline size_t extract_node_id(const std::string& s, size_t length_s, size_t& i) {
    size_t node_id = 0;
    while (i < length_s && s[i] != '<' && s[i] != '>') {
        node_id = node_id * 10 + (s[i] - '0');
        i++;
    }
    return node_id;
}

// Decompose path Path_traversal_t to vector Edge_t
std::vector<stoat_vcf::Edge_t> decompose_path_to_edges(const stoat_vcf::Path_traversal_t& list_paths) {
    std::vector<stoat_vcf::Edge_t> edges;
    const std::vector<Node_traversal_t>& list_nodes = list_paths.get_paths();
    size_t length_s = list_nodes.size();
    edges.reserve(length_s - 1); // Reserve memory

    for (size_t i = 0; i < length_s - 1; ++i) {
        edges.emplace_back(list_nodes[i], list_nodes[i + 1]);
    }

    return edges;
}

// Decompose a list of paths Path_traversal_t into a vector of Edge_t
const std::vector<std::vector<stoat_vcf::Edge_t>> decompose_path_list_path(const std::vector<stoat_vcf::Path_traversal_t>& list_paths) {
    size_t size_list_paths = list_paths.size();
    std::vector<std::vector<stoat_vcf::Edge_t>> paths_snarl;
    for (const stoat_vcf::Path_traversal_t& path : list_paths) {
        paths_snarl.push_back(decompose_path_to_edges(path));
    }
    return paths_snarl;
}

// Decompose path std::string to vector Edge_t
std::vector<stoat_vcf::Edge_t> decompose_path_str_to_edge(const std::string& s) {
    std::vector<stoat_vcf::Edge_t> edges;
    std::vector<Node_traversal_t> nodes;

    size_t i = 0;
    while (i < s.size()) {
        if (s[i] == '>' || s[i] == '<') {
            bool is_rev = (s[i] == '<');
            ++i;

            size_t node_id = 0;
            while (i < s.size() && isdigit(s[i])) {
                node_id = node_id * 10 + (s[i] - '0');
                ++i;
            }
            nodes.emplace_back(node_id, is_rev);
        } else {
            ++i; // Skip invalid characters
        }
    }

    for (size_t j = 0; j + 1 < nodes.size(); ++j) {
        edges.emplace_back(nodes[j], nodes[j + 1]);
    }

    return edges;
}

// Decompose a list of paths str into a vector of Edge_t
const std::vector<std::vector<stoat_vcf::Edge_t>> decompose_path_list_str(const std::vector<std::string>& list_paths) {
    std::vector<std::vector<stoat_vcf::Edge_t>> paths_snarl;
    for (const auto& path : list_paths) {
        paths_snarl.push_back(decompose_path_str_to_edge(path));
    }
    return paths_snarl;
}

// Function to identify the path in the edge matrix
std::vector<size_t> identify_path(
    const std::vector<Edge_t>& list_edge_path,
    const stoat_vcf::EdgeBySampleMatrix& edge_matrix,
    const size_t num_cols) {

    std::vector<size_t> rows_to_check;
    rows_to_check.reserve(list_edge_path.size());

    // Map snarl names to row indices
    for (const Edge_t& edge : list_edge_path) {
        const auto& [node_id_1, node_id_2] = edge.print_pair_edge(); // Convert Edge_t to std::pair<size_t, size_t>
        
        // Skip if snarl contains '*' (here * == 0) aka complex path
        if (node_id_1 == 0 || node_id_2 == 0) {
            continue;
        }
        size_t row_index = edge_matrix.find_edge(edge);
        if (row_index != std::numeric_limits<size_t>::max()) {
            rows_to_check.push_back(row_index);
        } else {
            return {}; // If any snarl isn't found, abort early
        }
    }

    std::vector<size_t> idx_srr_save;
    idx_srr_save.reserve(num_cols);

    // Loop columns first (better cache locality if matrix is column-major or similar)
    for (size_t col = 0; col < num_cols; ++col) {
        bool all_ones = true;
        for (size_t row : rows_to_check) {
            if (!edge_matrix(row, col)) {
                all_ones = false;
                break;
            }
        }
        if (all_ones) {
            idx_srr_save.push_back(static_cast<int>(col));
        }
    }
    return idx_srr_save;
}


void BinarySnarlAnalyzer::analyze_and_write_snarl(const std::string& chr, 
                                                  const Snarl_data_t& snarl_data_s,
                                                  const std::string& regression_dir, 
                                                  std::ofstream& outf) {

        std::ostringstream oss;

        for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
            if (i != 0) oss << ",";
            oss << snarl_data_s.type_variants[i];
        }

        std::string type_var_str = oss.str();

        size_t length_column_headers = snarl_data_s.snarl_paths.size();
        std::vector<size_t> g0(length_column_headers, 0);
        std::vector<size_t> g1(length_column_headers, 0);

        size_t total_sum = stoat_vcf::create_binary_table(g0, g1, binary_phenotype, snarl_data_s.snarl_paths, length_column_headers, sample_count, edge_matrix);
        bool df_filtration = check_MAF_threshold_binary(g0, g1, total_sum, length_column_headers, maf);

        // Binary analysis single test
        if (!df_filtration) { // good df
            const auto& [group_paths, 
                allele_number_str, min_row_index_str, numb_colum_str, 
                inter_group_str, average_str] = stoat_vcf::binary_stat_test(g0, g1);

            const auto& [fastfisher_p_value, chi2_p_value] = fk.fisher_khi2(g0, g1);

            # pragma omp critical (outf) 
            {
                write_binary(outf, chr, snarl_data_s, type_var_str, fastfisher_p_value, chi2_p_value, "", allele_number_str, min_row_index_str,
                            numb_colum_str, inter_group_str, average_str, group_paths);
            }
        }
    }
}

void BinaryCovarSnarlAnalyzer::analyze_and_write_snarl(const std::string& chr, 
    const Snarl_data_t& snarl_data_s,
    const std::string& regression_dir, 
    std::ofstream& outf) {

    std::ostringstream oss;

    for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
        if (i != 0) oss << ",";
        oss << snarl_data_s.type_variants[i];
    }

    std::string type_var_str = oss.str();

    const auto& [df, phenotype_filtered, allele_number, allele_paths] = create_quantitative_table(sample_count, snarl_data_s.snarl_paths, binary_phenotype, edge_matrix);
    bool df_filtration = check_MAF_threshold_quantitative(df, maf);

    if (!df_filtration) { // filtred variant
        // logistic regression with covariates if not empty
        const auto& [p_value, beta, se, r2] = lr.logistic_regression(df, phenotype_filtered, covar);

        // Plot regression table
        if (table_threshold != -1 && stoat_vcf::isPValueSignificant(table_threshold, p_value)) {
            std::string variant_file_name = regression_dir + "/" + pairToString(snarl_data_s.snarl_ids) + ".tsv";
            stoat_vcf::writeSignificantTableToTSV(df,stringToVector<std::string>(vectorPathToString(snarl_data_s.snarl_paths)), edge_matrix.sampleNames, variant_file_name);
        }
        # pragma omp critical (outf) 
        {
            write_binary_covar(outf, chr, snarl_data_s, type_var_str, p_value, "", r2, beta, se, allele_number, allele_paths);
        }
    }
}

// Quantitative Table Generation
void QuantitativeSnarlAnalyzer::analyze_and_write_snarl(const std::string &chr, const Snarl_data_t& snarl_data_s,
    const std::string& regression_dir, std::ofstream& outf) {

    const auto& [df, phenotype_filtered, allele_number, allele_paths] = create_quantitative_table(sample_count, snarl_data_s.snarl_paths, quantitative_phenotype, edge_matrix);
    bool df_filtration = check_MAF_threshold_quantitative(df, maf);
    
    // make a std::string separated by ',' from a vector of std::string
    std::ostringstream oss;
    for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
        if (i != 0) oss << ","; // Add comma before all elements except the first
        oss << snarl_data_s.type_variants[i];
    }

    std::string type_var_str = oss.str();
    std::stringstream data;
    
    if (!df_filtration) { // filtred variant
        auto [p_value, beta, se, r2] = lr.linear_regression(df, phenotype_filtered, covar);
        
        if (table_threshold != -1 && stoat_vcf::isPValueSignificant(table_threshold, p_value)) {
            std::string variant_file_name = regression_dir + "/" + pairToString(snarl_data_s.snarl_ids) + ".tsv";
            stoat_vcf::writeSignificantTableToTSV(df,stringToVector<std::string>(vectorPathToString(snarl_data_s.snarl_paths)), edge_matrix.sampleNames, variant_file_name);
        }
        
        #pragma omp critical (outf)
        {
            write_quantitative(outf, chr, snarl_data_s, type_var_str, p_value, "", r2, beta, se, allele_number, allele_paths);
        }
    }
}

// Identify genes index that will be tested for this snarl by matching position
// eqtl : <gene_name, gene_expression, start_pos, end_pos>
std::vector<size_t> found_gene_snarl(
    const std::vector<Qtl_data>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold) {

    std::vector<size_t> gene_index;
    size_t start_pos_threshold = (start_pos > windows_gene_threshold) ? start_pos - windows_gene_threshold : 0;
    size_t end_pos_threshold = end_pos + windows_gene_threshold;

    for (size_t i = 0; i < gene_position.size(); ++i) {
        size_t gene_start = gene_position[i].start_pos;
        size_t gene_end = gene_position[i].end_pos;

        // Check if the gene overlaps with the snarl region
        if (!(gene_end < start_pos_threshold || gene_start > end_pos_threshold)) {
            gene_index.push_back(i);
        }
    }
    return gene_index;
}

void EQTLSnarlAnalyzer::analyze_and_write_snarl(const Snarl_data_t& snarl_data_s) {

    std::vector<size_t> list_gene_index = found_gene_snarl(eqtl, snarl_data_s.start_positions, snarl_data_s.end_positions, windows_gene_threshold);
    const auto& [df, index_filtered, allele_number, allele_paths] = create_eqtl_table(sample_count, snarl_data_s.snarl_paths, edge_matrix);
    bool df_filtration = check_MAF_threshold_quantitative(df, maf);

    for (size_t i = 0; i < list_gene_index.size(); ++i) {
        size_t gene_idx = list_gene_index[i];
        std::string gene_name = eqtl[gene_idx].geneName;
        std::vector<double> gene_expression = eqtl[gene_idx].sampleExpresion;
        retain_indices(gene_expression, index_filtered);

        // make a std::string separated by ',' from a vector of std::string
        std::ostringstream oss;
        for (size_t i = 0; i < snarl_data_s.type_variants.size(); ++i) {
            if (i != 0) oss << ","; // Add comma before all elements except the first
            oss << snarl_data_s.type_variants[i];
        }

        std::string type_var_str = oss.str();
        std::stringstream data;

        if (!df_filtration) { // filtred variant
            auto [p_value, beta, se, r2] = lr.linear_regression(df, gene_expression, covar);

            if (table_threshold != -1 && stoat_vcf::isPValueSignificant(table_threshold, p_value)) {
                std::string variant_file_name = regression_dir + "/" + pairToString(snarl_data_s.snarl_ids) + ".tsv";
                stoat_vcf::writeSignificantTableToTSV(df,stringToVector<std::string>(vectorPathToString(snarl_data_s.snarl_paths)), edge_matrix.sampleNames, variant_file_name);
            }

            #pragma omp critical (outf)

            {
                stoat_vcf::write_eqtl(outf, chr, snarl_data_s, type_var_str, gene_name, p_value, "", r2, beta, se, allele_number, allele_paths);
            }
        }
    }
}

bool check_MAF_threshold_quantitative(const std::vector<std::vector<double>>& df, const double& maf) {    
    
    if ((df.size() < 2)) {return true;}

    double totalSum = 0;
    size_t numPaths = df[0].size(); // Get the number of paths from the first element
    std::vector<double> table(numPaths, 0); // Initialize vector with the correct size

    // Compute total sum of all elements in the matrix
    for (const auto& vector : df) {
        for (size_t i = 0; i < vector.size(); i++) {
            table[i] += vector[i];
            totalSum += vector[i];
        }
    }
}

bool check_MAF_threshold_binary(
    const std::vector<size_t>& g0, 
    const std::vector<size_t>& g1,
    size_t totalSum, 
    size_t length_column_headers, 
    double maf) {

    if (g0.empty() || g1.empty()) return true;

    for (size_t i = 0; i < length_column_headers; ++i) {
        size_t columnSum = g0[i] + g1[i];
        if (static_cast<double>(columnSum) / totalSum >= maf) {
            return true;  // MAF threshold met
        }
    }

    return false;  // No column met MAF threshold
}

} // end namespace stoat_vcf
