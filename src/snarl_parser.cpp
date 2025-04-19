#include "snarl_parser.hpp"
#include "matrix.hpp"
#include "binary_analysis.hpp"
#include "quantitative_analysis.hpp"
#include "utils.hpp"
#include "lmm.hpp"
#include "arg_parser.hpp"

using namespace std;

void chromosome_chuck_make_bed(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    string& output_dir) {

    const std::string output_bed = output_dir + ".bed";
    const std::string output_bim = output_dir + ".bim";

    std::ofstream outbim(output_bim);
    std::ofstream outbed(output_bed, std::ios::binary);  // Open BED file as binary
    
    // Write the 3-byte 'BED' header for the BED file
    char bed_magic[] = {0x6C, 0x1B, 0x01};  // PLINK header: 0x6C ('l'), 0x1B, 0x01 (snp-major mode)
    outbed.write(bed_magic, 3);
    
    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << "> " << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr[chr];

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

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    const unordered_map<string, double>& pheno, std::unordered_map<std::string, std::vector<double>> covar,
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const std::string& output_quantitive) {

    std::ofstream outf(output_quantitive, std::ios::binary);
    std::string headers;
    if (covar.size() > 0) {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP\tRSQUARED\tBETA\tSE\tALLELE_NUM\n";
    } else {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tRSQUARE\tBETA\tSE\tALLELE_NUM\n";
    }
    outf.write(headers.c_str(), headers.size());

    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << "> " << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.quantitative_table(snarl, pheno, chr, covar, maf, kinship, num_threads, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

// void chromosome_chuck_eqtl(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
//     const std::vector<std::string> &list_samples,
//     unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
//     const vector<QTL> qtl_pheno, std::ofstream& outf) {

//     std::cout << "GWAS analysis for chromosome : " << std::endl;
//     while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

//         string chr = bcf_hdr_id2name(hdr, rec->rid);
//         std::cout << "> " << chr << std::endl;
//         size_t size_chr = snarl_chr[chr].size();

//         // Make genotype matrix by chromosome    
//         auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
//         ptr_vcf = ptr_vcf_new;
//         hdr = hdr_new;
//         rec = rec_new;

//         auto& snarl = snarl_chr[chr];

//         // Gwas analysis by chromosome
//         // vcf_object.eqtl_table(snarl, pheno, chr, outf);
//     }
//     // Cleanup
//     bcf_destroy(rec);
//     bcf_hdr_destroy(hdr);
//     bcf_close(ptr_vcf);
// }

void chromosome_chuck_binary(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples, 
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    const unordered_map<string, bool>& pheno, std::unordered_map<std::string, std::vector<double>> covar, 
    const double& maf, const KinshipMatrix& kinship, 
    const size_t& num_threads, const std::string& output_binary) {

    std::ofstream outf(output_binary, std::ios::binary);
    std::string headers;
    if (covar.size() > 0) {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP\tP_ADJUSTED\tBETA\tSE\tALLELE_NUM\n";
    } else {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tP_ADJUSTED\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n";
    }
    outf.write(headers.c_str(), headers.size());

    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << "> " << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;
        auto& snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.binary_table(snarl, pheno, chr, covar, maf, kinship, num_threads, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

SnarlParser::SnarlParser(const vector<string>& sample_names, size_t num_paths_chr) : 
    sampleNames(sample_names), matrix(num_paths_chr*4, sample_names.size() * 2)
{}

std::pair<std::vector<size_t>, std::vector<size_t>> SnarlParser::create_table_short_path(const vector<std::string>& list_path_snarl) {

    size_t length_column = list_path_snarl.size();
    std::vector<size_t> allele_number_list(length_column, 0);
    size_t length_sample = sampleNames.size(); // get from the SnarlParser object

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<size_t>> genotypes(length_sample, std::vector<size_t>(length_column, 0));

    // Genotype paths
    for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
        const std::string& path_snarl = list_path_snarl[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, matrix, length_sample*2);

        for (auto idx : idx_srr_save) {
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

    return {genotypes_transposed[major_index_1], genotypes_transposed[major_index_1]};
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

void SnarlParser::create_bim_bed(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls, 
                                string chromosome, std::ofstream& outbim, std::ofstream& outbed) {

    // Iterate over each snarl
    // <snarl, paths, pos, type>
    for (const auto& [snarl, list_snarl, position, type] : snarls) {

        // if (list_snarl.size() > 2) {continue;} // avoid multiallelic var
        const size_t sample_number = sampleNames.size();  // Number of individuals

        // Generate a genotype table for this snarl
        auto [allele_vector_0, allele_vector_1] = create_table_short_path(list_snarl);
        
        std::string allele1 = "A";  // Placeholder for allele 1
        std::string allele2 = "T";  // Placeholder for allele 2

        // chr id genetic_distance pos allele1 allele2
        outbim << chromosome << "\t" << snarl << "\t0\t" << position
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

// Function to extract an integer from a string starting at index `i`
std::pair<int, std::string> determine_str(const std::string& s, size_t length_s, size_t i) {
    int start_idx = i;
    while (i < length_s && s[i] != '>' && s[i] != '<') {
        i++;
    }
    return {i, s.substr(start_idx, i - start_idx)};
}

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s) {
    std::vector<std::string> result;
    int i = 0;
    int length_s = s.length();
    std::string prev_int, prev_sym;

    while (i < length_s) {
        char start_sym = s[i];
        i++;
        auto [new_i, current_int] = determine_str(s, length_s, i);
        i = new_i;

        if (!prev_int.empty() && !prev_sym.empty()) {
            result.push_back(prev_sym + prev_int + start_sym + current_int);
        }

        prev_int = current_int;
        prev_sym = start_sym;
    }
    return result;
}

// Function to decompose a list of snarl strings
const std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst) {
    std::vector<std::vector<std::string>> decomposed_list;
    for (const auto& s : lst) {
        decomposed_list.push_back(decompose_string(s));
    }
    return decomposed_list;
}

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<std::string, size_t>& orderedMap, const std::string& key, size_t lengthOrderedMap) {
    auto it = orderedMap.find(key);
    if (it != orderedMap.end()) {
        return it->second;
    } else {
        size_t newIndex = lengthOrderedMap;
        orderedMap[key] = newIndex;
        return newIndex;
    }
}

// Add True to the matrix if snarl is found
void SnarlParser::push_matrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn) {
    
    size_t lengthOrderedMap = rowHeaderDict.size();
    size_t idxSnarl = getOrAddIndex(rowHeaderDict, decomposedSnarl, lengthOrderedMap);
    size_t currentRowsNumber = matrix.getMaxElement();
    
    if (lengthOrderedMap > currentRowsNumber - 1) {
        matrix.expandMatrix();
    }

    matrix.set(idxSnarl, indexColumn);
}

// Function to parse VCF and fill matrix genotypes
std::tuple<SnarlParser, htsFile*, bcf_hdr_t*, bcf1_t*> make_matrix(htsFile *ptr_vcf, bcf_hdr_t *hdr, bcf1_t *rec, const std::vector<std::string> &sampleNames, string &chr, size_t &num_paths_chr) {

    SnarlParser snarl_parser(sampleNames, num_paths_chr);
    std::unordered_map<std::string, size_t> row_header_dict;

    // loop over the VCF file for each line and stop where chr is different
    while ((bcf_read(ptr_vcf, hdr, rec) >= 0) && (chr == bcf_hdr_id2name(hdr, rec->rid))) {
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

        std::string at_value(at_str);  // Convert C-string to C++ string
        free(at_str);  // Free HTSlib-allocated memory

        // Split by comma
        std::stringstream ss(at_value);
        std::string item;
        while (std::getline(ss, item, ',')) {
            path_list.push_back(item);
        }

        // Decompose snarl paths
        const std::vector<std::vector<std::string>> list_list_decomposed_snarl = decompose_snarl(path_list);
        
        for (int i = 0; i < rec->n_sample; ++i) {
            int allele_1 = bcf_gt_allele(gt[i * 2]);
            int allele_2 = bcf_gt_allele(gt[i * 2 + 1]);
            size_t col_idx = i * 2;

            if (allele_1 != -1) { // Handle non-missing genotypes
                for (const auto &decompose_allele_1 : list_list_decomposed_snarl[allele_1]) {
                    snarl_parser.push_matrix(decompose_allele_1, row_header_dict, col_idx);
                }
            }
            if (allele_2 != -1) { // Handle non-missing genotypes
                for (const auto &decompose_allele_2 : list_list_decomposed_snarl[allele_2]) {
                    snarl_parser.push_matrix(decompose_allele_2, row_header_dict, col_idx + 1);
                }
            }
        }

        free(gt);
    }

    snarl_parser.matrix.set_row_header(row_header_dict);
    snarl_parser.matrix.shrink(row_header_dict.size());
    snarl_parser.matrix.set_end_dict();
    return std::make_tuple(snarl_parser, ptr_vcf, hdr, rec);
}

std::vector<int> identify_correct_path(
    const std::vector<std::string>& decomposed_snarl,
    const Matrix& matrix,
    const size_t num_cols) {

    std::vector<size_t> rows_to_check;
    rows_to_check.reserve(decomposed_snarl.size());

    // Map snarl names to row indices
    for (const auto& snarl : decomposed_snarl) {
        if (snarl.find("*") != std::string::npos) {
            continue;
        }
        auto it = matrix.find_snarl(snarl);
        if (it != matrix.get_end_dict()) {
            rows_to_check.push_back(it->second);
        } else {
            return {}; // If any snarl isn't found, abort early
        }
    }

    std::vector<int> idx_srr_save;
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

void SnarlParser::binary_table(const std::vector<std::tuple<std::string, std::vector<std::string>, std::string, std::vector<std::string>>>& snarls,
                               const std::unordered_map<std::string, bool>& binary_phenotype, const std::string& chr,
                               const std::unordered_map<std::string, std::vector<double>>& covar,
                               const double& maf, const KinshipMatrix& kinship, 
                               const size_t& num_threads, std::ofstream& outf) {

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

            for (size_t itr = start; itr < end; ++itr) {
                const auto& [snarl, list_snarl, pos, type_var] = snarls[itr];
                
                std::ostringstream oss;
                for (size_t i = 0; i < type_var.size(); ++i) {
                    if (i != 0) oss << ",";
                    oss << type_var[i];
                }

                std::string type_var_str = oss.str();
                std::stringstream data;

                if (!covar.empty()) {
                    // Logistic regression
                    auto [df, allele_number] = create_quantitative_table(sampleNames, list_snarl, matrix);
                    bool df_filtration = false;
                    bool df_empty = false;
            
                    if (allele_number < 2) {
                        df_empty = true;
                    } else {
                        df_filtration = check_MAF_threshold_quantitative(df, maf);
                    }
                    
                    // chr, pos, snarl, type, p_value, p_adjusted, r2, beta, se, allele_number
                    if (df_empty || !df_filtration) {
                        data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
                        << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA"
                        << "\t" << "NA" << "\t" << allele_number << "\n";
                    } else {
                        const auto& [beta, se, p_value] = logistic_regression(df, binary_phenotype, covar);
            
                        // chr, pos, snarl, type, p_value, p_adjusted, t-dist, beta, se, allele_number
                        data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
                        << "\t" << p_value << "\t" << "" << "\t" << beta << "\t" << se 
                        << "\t" << allele_number << "\n";
                    }

                } else {
                    size_t length_column_headers = list_snarl.size();
                    std::vector<size_t> g0(length_column_headers, 0); // can be replace by size_t arr[length_column_headers] = {0};
                    std::vector<size_t> g1(length_column_headers, 0); // can be replace by size_t arr[length_column_headers] = {0};
                    bool df_filtration = create_binary_table(g0, g1, binary_phenotype, list_snarl, sampleNames, matrix, maf);

                    std::string fastfisher_p_value = "NA", chi2_p_value = "NA",
                    group_paths = "NA", allele_number_str = "NA", min_row_index_str = "NA",
                    numb_colum_str = "NA", inter_group_str = "NA", average_str = "NA";
    
                    // Binary analysis single test
                    if (!df_filtration) { // good df
                        binary_stat_test(g0, g1, fastfisher_p_value, chi2_p_value, group_paths,
                            allele_number_str, min_row_index_str, numb_colum_str, inter_group_str, average_str);
                    }

                    data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
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
void SnarlParser::quantitative_table(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls,
                                        const std::unordered_map<std::string, double>& quantitative_phenotype, const string &chr,
                                        const std::unordered_map<std::string, std::vector<double>>& covar,
                                        const double& maf, const KinshipMatrix& kinship, const size_t& num_threads, std::ofstream& outf) {

    // Iterate over each snarl
    for (size_t itr = 0; itr < snarls.size(); ++itr) {
        const auto& [snarl, list_snarl, pos, type_var] = snarls[itr];

        auto [df, allele_number] = create_quantitative_table(sampleNames, list_snarl, matrix);
        bool df_filtration = false;
        bool df_empty = false;

        if (allele_number < 2) {
            df_empty = true;
        } else {
            df_filtration = check_MAF_threshold_quantitative(df, maf);
        }

        // make a string separated by ',' from a vector of string
        std::ostringstream oss;
        for (size_t i = 0; i < type_var.size(); ++i) {
            if (i != 0) oss << ","; // Add comma before all elements except the first
            oss << type_var[i];
        }
        std::string type_var_str = oss.str();
        std::stringstream data;

        // chr, pos, snarl, type, p_value, p_adjusted, r2, beta, se, allele_number
        if (df_empty || !df_filtration) {
            data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
            << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA"
            << "\t" << "NA" << "\t" << allele_number << "\n";
         
        } else if (covar.size() > 0) { // lmm
            const auto& [t_dist, beta, se, p_value] = lmm_quantitative(df, quantitative_phenotype, kinship, covar);

            // chr, pos, snarl, type, p_value, p_adjusted, t-dist, beta, se, allele_number
            data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
            << "\t" << p_value << "\t" << "" << "\t" << t_dist << "\t" 
            << beta << "\t" << se << "\t" << allele_number << "\n";

        } else { // single test
            const auto& [r2, beta, se, p_value] = linear_regression(df, quantitative_phenotype);

            // chr, pos, snarl, type, p_value, p_adjusted, r2, beta, se, allele_number
            data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
            << "\t" << p_value  << "\t" << "" << "\t" <<r2 << "\t" << beta << "\t" << se 
            << "\t" << allele_number << "\n";
        }

        outf.write(data.str().c_str(), data.str().size());
    }
}

bool check_MAF_threshold_quantitative(const std::unordered_map<std::string, std::vector<size_t>>& df, const double& maf) {    
    int totalSum = 0;
    size_t numPaths = df.begin()->second.size(); // Get the number of paths from the first element
    std::vector<int> table(numPaths, 0); // Initialize vector with the correct size

    // Compute total sum of all elements in the matrix
    for (const auto& [key, paths] : df) {
        for (size_t i = 0; i < paths.size(); i++) {
            table[i] += paths[i];
            totalSum += paths[i];
        }
    }
    
    // Check if any column's sum proportion exceeds the threshold
    for (int val : table) {
        if (static_cast<double>(val) / totalSum >= maf) {
            return false; // If any value exceeds the threshold, return false
        }
    }
    
    return true; // If all values are within the threshold, return true
}

std::unordered_map<std::string, std::vector<double>> convertBinaryGroups(
    const std::unordered_map<std::string, bool>& binary_phenotype) {

    std::unordered_map<std::string, std::vector<double>> converted_map;

    for (const auto& entry : binary_phenotype) {
        const std::string& sample_id = entry.first;
        bool group_value = entry.second;

        // Convert the binary value (bool) to a vector of double (1.0 or 0.0)
        std::vector<double> group_vector = {group_value ? 1.0 : 0.0};

        // Store the converted vector in the result map
        converted_map[sample_id] = group_vector;
    }

    return converted_map;
}
