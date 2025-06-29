#include "snarl_analyser.hpp"
#include "matrix.hpp"
#include "binary_test.hpp"
#include "quatitative_test.hpp"
#include "utils.hpp"
#include "lmm.hpp"
#include "arg_parser.hpp"

namespace stoat_vcf {

void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples, 
    const unordered_map<string, std::vector<Snarl_data_t>> &snarl_chr,
    const std::vector<bool>& binary_pheno, std::vector<std::vector<double>> covar, 
    const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
    const double& table_threshold, const std::string& regression_dir,
    const std::string& output_binary) {

    std::ofstream outf(output_binary, std::ios::binary);
    std::string headers;
    if (covar.size() > 0) {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS\n";
    } else {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tP_ADJUSTED\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n";
    }
    outf.write(headers.c_str(), headers.size());

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
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;
        auto& snarl = snarl_chr.at(chr);

        cout << "make_matrix done" << std::endl;
        // Gwas analysis by chromosome
        vcf_object.binary_table(snarl, binary_pheno, chr, covar, maf, kinship, num_threads, table_threshold, regression_dir, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}


void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const unordered_map<string, std::vector<Snarl_data_t>> &snarl_chr,
    const std::vector<double>& quantitative_phenotype, std::vector<std::vector<double>> covar,
    const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
    const double& table_threshold, const std::string& regression_dir,
    const std::string& output_quantitive) {

    std::ofstream outf(output_quantitive, std::ios::binary);
    std::string headers;
    headers = "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS\n";
    outf.write(headers.c_str(), headers.size());

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
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr.at(chr);

        // Gwas analysis by chromosome
        vcf_object.quantitative_table(snarl, quantitative_phenotype, chr, covar, maf, kinship, num_threads, table_threshold, regression_dir, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}


void chromosome_chuck_eqtl(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    const std::unordered_map<std::string, std::vector<Snarl_data_t>> &snarl_chr,
    const std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>>& eqtl_map,
    const std::vector<std::vector<double>>& covar,
    const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
    const double& table_threshold, const std::string& regression_dir,
    const size_t& windows_gene_threshold, const std::string& out_eqtl) {

    std::ofstream outf(out_eqtl, std::ios::binary);
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tGENE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\tALLELE_PATHS\n";
    outf.write(headers.c_str(), headers.size());
    
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
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr.at(chr);
        auto& eqtl = eqtl_map.at(chr);

        // Gwas analysis by chromosome
        vcf_object.eqtl_table(snarl, eqtl, chr, covar, maf, kinship, num_threads, table_threshold, regression_dir, windows_gene_threshold, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
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
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr.at(chr);

        // Gwas analysis by chromosome
        vcf_object.create_bim_bed(snarl, chr, outbim, outbed);
    }
    
    // Cleanup
    outbim.close();
    outbed.close();
    
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

SnarlAnalyser::SnarlAnalyser(const std::vector<std::string>& sample_names, size_t num_paths_chr) : 
    sampleNames(sample_names), matrix(num_paths_chr*4, sample_names.size() * 2)
{}

std::pair<std::vector<size_t>, std::vector<size_t>> SnarlAnalyser::create_table_short_path(const std::vector<Path_traversal_t>& list_path_snarl) {
    // TODO: Could just make the transposed version to start with?

    size_t length_column = list_path_snarl.size();
    std::vector<size_t> allele_number_list(length_column, 0);
    size_t length_sample = sampleNames.size(); // get from the SnarlAnalyser object

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<size_t>> genotypes(length_sample, std::vector<size_t>(length_column, 0));

    // Genotype paths
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const Path_traversal_t& path_snarl = list_path_snarl[col_idx];
        std::vector<Edge_t> decomposed_snarl = decompose_path_to_edges(path_snarl);

        // Identify correct paths
        std::vector<size_t> idx_srr_save = identify_path(decomposed_snarl, matrix, length_sample*2);

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

void SnarlAnalyser::create_bim_bed(const std::vector<Snarl_data_t>& snarls, 
                                std::string chromosome, std::ofstream& outbim, std::ofstream& outbed) {

    // Iterate over each snarl
    // <snarl, paths, pos, type>
    for (const Snarl_data_t& snarl_data_s : snarls) {

        std::string snarl_id = pairToString(snarl_data_s.get_snarl_id());
        std::vector<Path_traversal_t> list_path_snarl = snarl_data_s.get_snarl_paths();
        size_t start_pos = snarl_data_s.get_start_positions();

        // if (list_path_snarl.size() > 2) {continue;} // avoid multiallelic var
        const size_t sample_number = sampleNames.size();  // Number of individuals

        // Generate a genotype table for this snarl
        auto [allele_vector_0, allele_vector_1] = create_table_short_path(list_path_snarl);
        
        std::string allele1 = "A";  // Placeholder for allele 1
        std::string allele2 = "T";  // Placeholder for allele 2

        // chr id genetic_distance pos allele1 allele2
        outbim << chromosome << "\t" << snarl_id << "\t0\t" << start_pos
                << "\t" << allele1 << "\t" << allele2 << "\n";
        
        // Write the genotypes for this SNP to the BED file
        unsigned char packed_byte = 0;  // A byte to store genotypes of 4 individuals
        int bit_pos = 0;

        // Loop through each sample (pair of alleles per individual)
        for (size_t snarl_list_idx = 0; snarl_list_idx < sample_number; ++snarl_list_idx) {

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
std::vector<Edge_t> decompose_path_to_edges(const Path_traversal_t& list_paths) {
    std::vector<Edge_t> edges;
    const std::vector<Node_traversal_t>& list_nodes = list_paths.get_paths();
    size_t length_s = list_nodes.size();
    edges.reserve(length_s - 1); // Reserve memory

    for (size_t i = 0; i < length_s - 1; ++i) {
        edges.emplace_back(list_nodes[i], list_nodes[i + 1]);
    }

    return edges;
}

// Decompose a list of paths Path_traversal_t into a vector of Edge_t
const std::vector<std::vector<Edge_t>> decompose_path_list_path(const std::vector<Path_traversal_t>& list_paths) {
    size_t size_list_paths = list_paths.size();
    std::vector<std::vector<Edge_t>> paths_snarl;
    for (const Path_traversal_t& path : list_paths) {
        paths_snarl.push_back(decompose_path_to_edges(path));
    }
    return paths_snarl;
}

// Decompose path std::string to vector Edge_t
std::vector<Edge_t> decompose_path_str_to_edge(const std::string& s) {
    std::vector<Edge_t> edges;
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
const std::vector<std::vector<Edge_t>> decompose_path_list_str(const std::vector<std::string>& list_paths) {
    std::vector<std::vector<Edge_t>> paths_snarl;
    for (const auto& path : list_paths) {
        paths_snarl.push_back(decompose_path_str_to_edge(path));
    }
    return paths_snarl;
}

// Retrieve the index of `key` if it exists in the dict. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<Edge_t, size_t>& edge_index_dict, const Edge_t& key, const size_t& size_edge_index_dict) {
    auto it = edge_index_dict.find(key);
    if (it != edge_index_dict.end()) {
        return it->second;
    } else {
        size_t newIndex = size_edge_index_dict;
        edge_index_dict[key] = newIndex;
        return newIndex;
    }
}

// Add True to the matrix if snarl is found
void SnarlAnalyser::push_matrix(const Edge_t& EdgePath, std::unordered_map<Edge_t, size_t>& edge_index_dict, size_t indexColumn) {
    
    size_t lengthOrderedMap = edge_index_dict.size();
    size_t idxSnarl = getOrAddIndex(edge_index_dict, EdgePath, lengthOrderedMap);
    size_t currentRowsNumber = matrix.getMaxElement();

    if (lengthOrderedMap > currentRowsNumber - 1) {
        matrix.expandMatrix();
    }

    matrix.set(idxSnarl, indexColumn);
}

// Function to parse VCF and fill matrix genotypes
std::tuple<SnarlAnalyser, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const std::vector<std::string> &sampleNames, std::string &chr, size_t &num_paths_chr) {

    SnarlAnalyser snarl_data(sampleNames, num_paths_chr);
    std::unordered_map<Edge_t, size_t> edge_dict;

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
        const std::vector<std::vector<Edge_t>> list_paths_edge = decompose_path_list_str(path_list);

        for (int i = 0; i < rec->n_sample; ++i) {
            int idex_path_allele_1 = bcf_gt_allele(gt[i * 2]);
            int idex_path_allele_2 = bcf_gt_allele(gt[i * 2 + 1]);
            size_t col_idx = i * 2;

            if (idex_path_allele_1 != -1) { // Handle missing genotypes
                for (const auto &edge_path_1 : list_paths_edge[idex_path_allele_1]) {
                    snarl_data.push_matrix(edge_path_1, edge_dict, col_idx);
                }
            }

            if (idex_path_allele_2 != -1) { // Handle missing genotypes
                for (const auto &edge_path_2 : list_paths_edge[idex_path_allele_2]) {
                    snarl_data.push_matrix(edge_path_2, edge_dict, col_idx + 1);
                }
            }
        }
        free(gt);

    } while ((bcf_read(ptr_vcf, hdr, rec) >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid)));

    snarl_data.matrix.set_row_header(edge_dict);
    snarl_data.matrix.shrink(edge_dict.size());
    snarl_data.matrix.set_end_dict();
    return std::make_tuple(snarl_data, ptr_vcf, hdr, rec);
}

// Function to identify the path in the edge matrix
std::vector<size_t> identify_path(
    const std::vector<Edge_t>& list_edge_path,
    const EdgeBySampleMatrix& matrix,
    const size_t num_cols) {

    std::vector<size_t> rows_to_check;
    rows_to_check.reserve(list_edge_path.size());

    // TODO: I think this is going through edges in the path through the snarl, not snarls
    // Map snarl names to row indices
    for (const Edge_t& edge : list_edge_path) {
        const auto& [node_id_1, node_id_2] = edge.print_pair_edge(); // Convert Edge_t to std::pair<size_t, size_t>
        
        // Skip if snarl contains '*' (here * == 0) aka complex path
        if (node_id_1 == 0 || node_id_2 == 0) {
            continue;
        }
        auto it = matrix.find_edge(edge);
        if (it != matrix.get_end_dict()) {
            rows_to_check.push_back(it->second);
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
            if (!matrix(row, col)) {
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

void SnarlAnalyser::binary_table(const std::vector<Snarl_data_t>& snarls,
                               const std::vector<bool>& binary_phenotype, const std::string& chr,
                               const std::vector<std::vector<double>>& covar,
                               const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
                               const double& table_threshold, const std::string& regression_dir, std::ofstream& outf) {

    size_t length_sample = sampleNames.size();
    const size_t total = snarls.size();
    size_t chunk_size = (total + num_threads - 1) / num_threads;
    std::mutex mutex_pvalues;
    std::mutex mutex_file;
    std::vector<std::thread> threads;

    cout << "Start binary_table " << std::endl;
    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        threads.emplace_back([&, thread_id]() {
            size_t start = thread_id * chunk_size;
            size_t end = std::min(start + chunk_size, total);
            std::stringstream local_buffer;

            for (size_t itr = start; itr < end; ++itr) {
                const Snarl_data_t& snarl_data_s = snarls[itr];
                cout << "start get snarl data" << std::endl;
                const auto& [snarl_id, list_path_snarl, start_pos, end_pos, type_var] = snarl_data_s.get_snarl();
                cout << "snarl_id : " << snarl_id << std::endl;

                std::ostringstream oss;
                for (size_t i = 0; i < type_var.size(); ++i) {
                    if (i != 0) oss << ",";
                    oss << type_var[i];
                }

                std::string type_var_str = oss.str();
                std::stringstream data;

                if (!covar.empty()) {
                    // Logistic regression
                    const auto& [df, phenotype_filtered, allele_number, allele_paths] = create_quantitative_table(length_sample, list_path_snarl, binary_phenotype, matrix);
                    bool df_filtration = false;
                    bool df_empty = false;

                    if (allele_number < 5) {
                        df_empty = true;
                    } else {
                        df_filtration = check_MAF_threshold_quantitative(df, maf);
                    }

                    std::string p_value = "", beta = "", se = "", r2 = "";

                    if (df_empty || df_filtration) { // filtred variant
                        // do not analyse this snarl
                        continue;
                    } else if (kinship.empty()) { // logistic regression + covar
                        logistic_regression(df, phenotype_filtered, p_value, beta, se, r2);
                        // glm_logistic_covar(df, phenotype_filtered, covar, p_value, beta, se, r2);
                    } else { // lmm
                        lmm_binary(df, phenotype_filtered, kinship, covar, p_value, beta, se, r2);
                    }
                    
                    // Plot regression table
                    if (table_threshold != -1 && stoat_vcf::isPValueSignificant(table_threshold, p_value)) {
                        std::string variant_file_name = regression_dir + "/" + snarl_id + ".tsv";
                        stoat_vcf::writeSignificantTableToTSV(df,stoat_vcf::stringToVector<std::string>(vectorPathToString(list_path_snarl)), sampleNames, variant_file_name);
                    }
    
                    // chr, pos, snarl, type, p_value, p_adjusted, t-dist, beta, se, allele_number
                    data << chr << "\t" << start_pos << "\t" << snarl_id << "\t" << type_var_str
                    << "\t" << p_value << "\t" << "" << "\t" << r2 << "\t" << beta << "\t" << se 
                    << "\t" << allele_number << "\t" << stoat_vcf::vectorToString(allele_paths) << "\n";
                                
                } else {
                    size_t length_column_headers = list_path_snarl.size();
                    size_t number_samples = sampleNames.size();
                    std::vector<size_t> g0(length_column_headers, 0); // can be replace by size_t arr[length_column_headers] = {0};
                    std::vector<size_t> g1(length_column_headers, 0); // can be replace by size_t arr[length_column_headers] = {0};

                    cout << "Creating binary table for snarl: " << snarl_id << std::endl;
                    size_t total_sum = create_binary_table(g0, g1, binary_phenotype, list_path_snarl, length_column_headers, number_samples, matrix);
                    bool df_filtration = check_MAF_threshold(g0, g1, total_sum, length_column_headers, maf);
                    cout << "Creating binary table finish" << std::endl;

                    std::string fastfisher_p_value = "NA", chi2_p_value = "NA",
                    group_paths = "NA", allele_number_str = "NA", min_row_index_str = "NA",
                    numb_colum_str = "NA", inter_group_str = "NA", average_str = "NA";

                    // Binary analysis single test
                    if (!df_filtration) { // good df
                        binary_stat_test(g0, g1, fastfisher_p_value, chi2_p_value, group_paths,
                            allele_number_str, min_row_index_str, numb_colum_str, inter_group_str, average_str);
                    }
                    
                    data << chr << "\t" << start_pos << "\t" << snarl_id << "\t" << type_var_str
                         << "\t" << fastfisher_p_value << "\t" << chi2_p_value << "\t" << ""
                         << "\t" << allele_number_str << "\t" << min_row_index_str << "\t" << numb_colum_str 
                         << "\t" << inter_group_str << "\t" << average_str << "\t" << group_paths << "\n";
                }

                local_buffer << data.str();
            }

            {
                std::lock_guard<std::mutex> lock(mutex_file);
                outf.write(local_buffer.str().c_str(), local_buffer.str().size());
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}

// Quantitative Table Generation
void SnarlAnalyser::quantitative_table(const std::vector<Snarl_data_t>& snarls,
                                       const std::vector<double>& quantitative_phenotype, 
                                       const std::string &chr,
                                       const std::vector<std::vector<double>>& covar,
                                       const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
                                       const double& table_threshold, const std::string& regression_dir, std::ofstream& outf) {

    size_t length_sample = sampleNames.size();
    const size_t total = snarls.size();
    size_t chunk_size = (total + num_threads - 1) / num_threads;
    std::mutex mutex_pvalues;
    std::mutex mutex_file;
    std::vector<std::thread> threads;

    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        threads.emplace_back([&, thread_id]() {
            size_t start = thread_id * chunk_size;
            size_t end = std::min(start + chunk_size, total);
            std::stringstream local_buffer;

            // Iterate over each snarl
            for (size_t itr = 0; itr < snarls.size(); ++itr) {
                const Snarl_data_t& snarl_data_s = snarls[itr];
                const auto& [snarl_id, list_path_snarl, start_pos, end_pos, type_var] = snarl_data_s.get_snarl();

                const auto& [df, phenotype_filtered, allele_number, allele_paths] = create_quantitative_table(length_sample, list_path_snarl, quantitative_phenotype, matrix);
                bool df_filtration = false;
                bool df_empty = false;

                if (allele_number < 5) {
                    df_empty = true;
                } else {
                    df_filtration = check_MAF_threshold_quantitative(df, maf);
                }

                // make a std::string separated by ',' from a vector of std::string
                std::ostringstream oss;
                for (size_t i = 0; i < type_var.size(); ++i) {
                    if (i != 0) oss << ","; // Add comma before all elements except the first
                    oss << type_var[i];
                }
                std::string type_var_str = oss.str();
                std::stringstream data;
                std::string p_value = "", beta = "", se = "", r2 = "";
                
                if (df_empty || df_filtration) { // filtred variant
                    // do not analyse this snarl
                    continue;
                } else if (covar.size() > 0 && !kinship.empty()) { // lmm
                    lmm_quantitative(df, phenotype_filtered, kinship, covar, p_value, beta, se, r2);

                } else if (covar.size() > 0 && kinship.empty()) { // glm
                    glm_quantitative(df, phenotype_filtered, covar, p_value, beta, se, r2);

                } else { // single test
                    linear_regression(df, phenotype_filtered, p_value, beta, se, r2);
                }
                
                if (table_threshold != -1 && stoat_vcf::isPValueSignificant(table_threshold, p_value)) {
                    std::string variant_file_name = regression_dir + "/" + snarl_id + ".tsv";
                    stoat_vcf::writeSignificantTableToTSV(df,stoat_vcf::stringToVector<std::string>(vectorPathToString(list_path_snarl)), sampleNames, variant_file_name);
                }
                
                // chr, pos, snarl, type, p_value, p_adjusted, r2, beta, se, allele_number
                data << chr << "\t" << start_pos << "\t" << snarl_id << "\t" << type_var_str
                << "\t" << p_value  << "\t" << "" << "\t" << r2 << "\t" << beta << "\t" << se 
                << "\t" << allele_number << "\t" << stoat_vcf::vectorToString(allele_paths) << "\n";
                local_buffer << data.str();
            }

            {
                std::lock_guard<std::mutex> lock(mutex_file);
                outf.write(local_buffer.str().c_str(), local_buffer.str().size());
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}

bool check_MAF_threshold_quantitative(const std::vector<std::vector<double>>& df, const double& maf) {    
    
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

    if (df[0].size() == 1) {
        if (totalSum == df.size()) {
            return true; // Case where all values are 1, return true avoiding const df 
        }
    
    } else {
        // Check if any column's sum proportion exceeds the threshold
        for (double val : table) {
            if (val / totalSum >= maf) {
                return true; // If any value exceeds the threshold, return false
            }
        }
    }

    return false; // If all values are within the threshold, return true
}

// Identify genes index that will be tested for this snarl by matching position
// eqtl : <gene_name, gene_expression, start_pos, end_pos>
std::vector<size_t> found_gene_snarl(
    const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& gene_position, 
    const size_t& start_pos, 
    const size_t& end_pos,
    const size_t& windows_gene_threshold) {

    std::vector<size_t> gene_index;
    size_t start_pos_threshold = (start_pos > windows_gene_threshold) ? start_pos - windows_gene_threshold : 0;
    size_t end_pos_threshold = end_pos + windows_gene_threshold;

    for (size_t i = 0; i < gene_position.size(); ++i) {
        size_t gene_start = std::get<2>(gene_position[i]);
        size_t gene_end = std::get<3>(gene_position[i]);

        // Check if the gene overlaps with the snarl region
        if (!(gene_end < start_pos_threshold || gene_start > end_pos_threshold)) {
            gene_index.push_back(i);
        }
    }
    return gene_index;
}

void SnarlAnalyser::eqtl_table(
    const std::vector<Snarl_data_t>& snarls,
    const std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>& eqtl,
    const std::string& chr, const std::vector<std::vector<double>>& covar,
    const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, 
    const double& table_threshold, const std::string& regression_dir, 
    const size_t& windows_gene_threshold, std::ofstream& outf) {

    size_t length_sample = sampleNames.size();
    const size_t total = snarls.size();
    size_t chunk_size = (total + num_threads - 1) / num_threads;
    std::mutex mutex_pvalues;
    std::mutex mutex_file;
    std::vector<std::thread> threads;

    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        threads.emplace_back([&, thread_id]() {
            size_t start = thread_id * chunk_size;
            size_t end = std::min(start + chunk_size, total);
            std::stringstream local_buffer;

            // Iterate over each snarl
            for (size_t itr = 0; itr < snarls.size(); ++itr) {
                const Snarl_data_t& snarl_data_s = snarls[itr];
                const auto& [snarl_id, list_path_snarl, start_pos, end_pos, type_var] = snarl_data_s.get_snarl();

                std::vector<size_t> list_gene_index = found_gene_snarl(eqtl, start_pos, end_pos, windows_gene_threshold);
                const auto& [df, index_filtered, allele_number, allele_paths] = create_eqtl_table(length_sample, list_path_snarl, matrix);
                bool df_filtration = false;
                bool df_empty = false;

                if (allele_number < 5) {
                    df_empty = true;
                } else {
                    if (df[0].size() > 1) {
                        df_filtration = check_MAF_threshold_quantitative(df, maf);
                    }
                }

                for (size_t i = 0; i < list_gene_index.size(); ++i) {
                    size_t gene_idx = list_gene_index[i];
                    std::string gene_name = std::get<0>(eqtl[gene_idx]);
                    std::vector<double> gene_expression = std::get<1>(eqtl[gene_idx]);
                    stoat_vcf::retain_indices(gene_expression, index_filtered);

                    // make a std::string separated by ',' from a vector of std::string
                    std::ostringstream oss;
                    for (size_t i = 0; i < type_var.size(); ++i) {
                        if (i != 0) oss << ","; // Add comma before all elements except the first
                        oss << type_var[i];
                    }

                    std::string type_var_str = oss.str();
                    std::stringstream data;
                    std::string p_value = "", beta = "", se = "", r2 = "";
    
                    if (df_empty || df_filtration) { // filtred variant
                        // do not analyse this snarl
                        continue;
                    } else if (covar.size() > 0 && !kinship.empty()) { // lmm
                        lmm_quantitative(df, gene_expression, kinship, covar, p_value, beta, se, r2);

                    } else if (covar.size() > 0 && kinship.empty()) { // glm
                        glm_quantitative(df, gene_expression, covar, p_value, beta, se, r2); // TODO : se nan problem

                    } else { // single test
                        linear_regression(df, gene_expression, p_value, beta, se, r2);
                    }

                    if (table_threshold != -1 && stoat_vcf::isPValueSignificant(table_threshold, p_value)) {
                        std::string variant_file_name = regression_dir + "/" + snarl_id + ".tsv";
                        stoat_vcf::writeSignificantTableToTSV(df,stoat_vcf::stringToVector<std::string>(vectorPathToString(list_path_snarl)), sampleNames, variant_file_name);
                    }

                   // "CHR\tPOS\tSNARL\tTYPE\tGENE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\n";
                    data << chr << "\t" << start_pos << "\t" << snarl_id << "\t" << type_var_str
                    << "\t" << gene_name << "\t" << p_value  << "\t" << "" << "\t" << r2
                    << "\t" << beta << "\t" << se << "\t" << allele_number << "\t" << stoat_vcf::vectorToString(allele_paths) << "\n";

                    local_buffer << data.str();
                }
            }

            {
                std::lock_guard<std::mutex> lock(mutex_file);
                outf.write(local_buffer.str().c_str(), local_buffer.str().size());
            }
        });
    }

    for (auto& t : threads) {
        t.join();
    }
}

} // end namespace stoat_vcf
