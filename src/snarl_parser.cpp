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
    const unordered_map<string, double>& pheno, string output_dir) {

    const std::string output_bed = output_dir + "genotype.bed";
    const std::string output_bim = output_dir + "genotype.bim";

    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.create_bim_bed(snarl, chr, output_bim, output_bed);
    }
    
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

void chromosome_chuck_quantitative(htsFile* &ptr_vcf, bcf_hdr_t* &hdr, bcf1_t* &rec, 
    const std::vector<std::string> &list_samples,
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> &snarl_chr,
    const unordered_map<string, double>& pheno, std::unordered_map<std::string, std::vector<double>> covar,
    const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
    const KinshipMatrix& kinship, const std::string& output_quantitive) {

    std::ofstream outf(output_quantitive, std::ios::binary);
    std::string headers;
    if (covar.size() > 0) {
        headers = "CHR\tPOS\tSNARL\tTYPE\tRSQUARED\tBETA\tSE\tP\tALLELE_NUM\n";
    } else {
        headers = "CHR\tPOS\tSNARL\tTYPE\tRSQUARE\tBETA\tSE\tP\tP_ADJUSTED\tALLELE_NUM\n";
    }
    outf.write(headers.c_str(), headers.size());

    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;

        auto& snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.quantitative_table(snarl, pheno, chr, covar, maf, pvalue_vector, kinship, outf);
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
//         std::cout << chr << std::endl;
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
    const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
    const KinshipMatrix& kinship, const std::string& output_binary) {

    std::ofstream outf(output_binary, std::ios::binary);
    std::string headers;
    if (covar.size() > 0) {
        headers = "CHR\tPOS\tSNARL\tTYPE\tBETA\tSE\tP\n";
    } else {
        headers = "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\tP_ADJUSTED\tALLELE_NUM\tMIN_ROW_INDEX\tNUM_COLUM\tINTER_GROUP\tAVERAGE\tGROUP_PATHS\n";
    }
    outf.write(headers.c_str(), headers.size());

    std::cout << "GWAS analysis for chromosome : " << std::endl;
    while (bcf_read(ptr_vcf, hdr, rec) >= 0) {

        string chr = bcf_hdr_id2name(hdr, rec->rid);
        std::cout << chr << std::endl;
        size_t size_chr = snarl_chr[chr].size();

        // Make genotype matrix by chromosome    
        auto [vcf_object, ptr_vcf_new, hdr_new, rec_new] = make_matrix(ptr_vcf, hdr, rec, list_samples, chr, size_chr);
        ptr_vcf = ptr_vcf_new;
        hdr = hdr_new;
        rec = rec_new;
        auto& snarl = snarl_chr[chr];

        // Gwas analysis by chromosome
        vcf_object.binary_table(snarl, pheno, chr, covar, maf, pvalue_vector, kinship, outf);
    }
    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(ptr_vcf);
}

SnarlParser::SnarlParser(const vector<string>& sample_names, size_t num_paths_chr) : 
    sampleNames(sample_names), matrix(num_paths_chr*4, sample_names.size() * 2)
{}

std::vector<int> SnarlParser::create_table_short_path(const std::string& list_path_snarl) {
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();

    // Iterate over each path_snarl in column_headers
    const std::string& path_snarl = list_path_snarl;
    const size_t number_sample = sampleNames.size();
    std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
    return identify_correct_path(decomposed_snarl, row_headers_dict, matrix, number_sample*2);
}

void SnarlParser::create_bim_bed(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls, 
                                string chromosome, const std::string& output_bim, const std::string& output_bed) {

    std::ofstream outbim(output_bim);
    std::ofstream outbed(output_bed, std::ios::binary);  // Open BED file as binary
    
    const size_t allele_number = sampleNames.size();  // Number of individuals

    if (!outbim.is_open() || !outbed.is_open()) {
        std::cerr << "Error opening output files!" << std::endl;
        return;
    }

    // Write the 3-byte 'BED' header for the BED file
    char bed_magic[] = {0x6C, 0x1B, 0x01};  // PLINK header: 0x6C ('l'), 0x1B, 0x01 (snp-major mode)
    outbed.write(bed_magic, 3);
    
    // Iterate over each snarl
    // <snarl, paths, pos, type>
    for (const auto& [snarl, list_snarl, position, type] : snarls) {

        if (list_snarl.size() > 2) {continue;} // avoid multiallelic var

        // Generate a genotype table for this snarl
        std::vector<int> table = create_table_short_path(list_snarl[0]);  // 2D vector, each row = SNP, each col = sample alleles

        std::string allele1 = "A";  // Placeholder for allele 1
        std::string allele2 = "T";  // Placeholder for allele 2

        // Write the BIM file line for the SNP
        outbim << chromosome << "\t" << snarl << "\t0\t" << position
                << "\t" << allele1 << "\t" << allele2 << "\n";
        
        // Write the genotypes for this SNP to the BED file
        unsigned char packed_byte = 0;  // A byte to store genotypes of 4 individuals
        int bit_pos = 0;

        // Loop through each sample (pair of alleles per individual)
        for (size_t snarl_list_idx = 0; snarl_list_idx < allele_number; ++snarl_list_idx) {
            int allele1 = table[2 * snarl_list_idx];      // First allele for the individual for the first paths
            int allele2 = table[2 * snarl_list_idx + 1];  // Second allele for the individual at SNP

            // Encode the genotype as a 2-bit value based on the alleles
            unsigned char encoded_genotype = 0;

            if (allele1 == allele2) {
                // Homozygous genotype (AA or aa)
                encoded_genotype = (allele1 == 0) ? 0b00 : 0b11;  // Example: 0 = AA, 1 = aa
            } else {
                // Heterozygous genotype (Aa)
                encoded_genotype = 0b01;
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
    outbim.close();
    outbed.close();
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
    while ((bcf_read(ptr_vcf, hdr, rec) >= 0) || (chr != bcf_hdr_id2name(hdr, rec->rid))) {
        bcf_unpack(rec, BCF_UN_STR);

        // Check the INFO field for LV (Level Variant) and skip if LV != 0
        int32_t *lv = nullptr;
        int n_lv = 0;
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

        if (nat > 0 && at_str) {
            std::string at_value(at_str);  // Convert C-string to C++ string
            free(at_str);  // Free HTSlib-allocated memory

            // Split by comma
            std::stringstream ss(at_value);
            std::string item;
            while (std::getline(ss, item, ',')) {
                path_list.push_back(item);
            }

        } else {
            std::cerr << "No AT field found" << std::endl;
            free(at_str);
        }

        // Decompose snarl paths
        const std::vector<std::vector<std::string>> list_list_decomposed_snarl = decompose_snarl(path_list);
        
        if (ngt > 0 && gt) {
            // loop over the genotypes
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
        } else {
            cerr << "No genotypes found" << std::endl;
        }
        free(gt);
    }

    snarl_parser.matrix.set_row_header(row_header_dict);
    snarl_parser.matrix.shrink(row_header_dict.size());
    return std::make_tuple(snarl_parser, ptr_vcf, hdr, rec);
}

std::vector<int> identify_correct_path(
    const std::vector<std::string>& decomposed_snarl,
    const std::unordered_map<std::string, size_t>& row_headers_dict,
    const Matrix& matrix,
    const size_t num_cols) {

    std::vector<int> rows_to_check;

    for (const auto& snarl : decomposed_snarl) {
        if (snarl.find("*") != std::string::npos) {
            continue; // Skip snarls containing "*"
        }
        auto it = row_headers_dict.find(snarl);
        if (it != row_headers_dict.end()) {
            rows_to_check.push_back(it->second);
        } else {
            return {}; // Return an empty vector if snarl is not in row_headers_dict
        }
    }

    // Check columns for all 1s in the specified rows
    std::vector<bool> columns_all_ones(num_cols, true);

    for (size_t col = 0; col < num_cols; ++col) {
        for (size_t row : rows_to_check) {
            if (!matrix(row, col)) {  // Use the `operator()` to access the matrix element
                columns_all_ones[col] = false;
                break; // Stop checking this column if any element is not 1
            }
        }
    }

    // Populate idx_srr_save with indices of columns where all elements are 1
    std::vector<int> idx_srr_save;
    for (size_t col = 0; col < num_cols; ++col) {
        if (columns_all_ones[col]) {
            idx_srr_save.push_back(col);
        }
    }
    return idx_srr_save;
}

// Binary Table Generation
void SnarlParser::binary_table(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls,
                               const std::unordered_map<std::string, bool>& binary_groups, const string &chr,
                               const std::unordered_map<std::string, std::vector<double>>& covar,
                               const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
                               const KinshipMatrix& kinship, std::ofstream& outf) {

    // Iterate over each snarl
    size_t itr = 0;
    for (const auto& tuple_snarl : snarls) {

        std::vector<std::string> list_snarl = std::get<1>(tuple_snarl);
        std::vector<std::vector<int>> df = create_binary_table(binary_groups, list_snarl, sampleNames, matrix);
        std::string snarl = std::get<0>(tuple_snarl), pos = std::get<2>(tuple_snarl);
        std::vector<std::string> type_var = std::get<3>(tuple_snarl);
        bool df_filtration = check_MAF_threshold_binary(df, maf);

        // make a string separated by , from a vector of string
        std::ostringstream oss;
            for (size_t i = 0; i < type_var.size(); ++i) {
                if (i != 0) oss << ","; // Add comma before all elements except the first
                oss << type_var[i];
            }

        std::string type_var_str = oss.str();
        std::vector<std::string> stats;
        std::stringstream data;
        
        if (covar.size() > 0) {
            // stats = lmm_binay(df, binary_phenotype, kinship, covar);
            // string p_value = df_filtration ? stats[2] : "NA";

            // // chr, pos, snarl, type variant, beta, se, p_value
            // data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
            // << "\t" << stats[0] << "\t" << stats[1] << "\t" << p_value << "\n";

        } else {
            // fastfisher_p_value, chi2_p_value, allele_number, min_row_index, numb_colum, inter_group, average, group_paths
            stats = binary_stat_test(df);
            string p_value_f, p_value_c; 

            if (df_filtration) {
                p_value_f = stats[0];
                p_value_c = stats[1];
            } else {
                p_value_f = p_value_c = "NA";
            }
            
            double mean_pvalue = mean_pvalue_from_strings(p_value_f, p_value_c);
            pvalue_vector.push_back({mean_pvalue, 1, itr});
            // chr, pos, snarl, type variant
            // fisher_p_value, chi2_p_value, adjusted_pvalue, allele_number, min_row_index, 
            // numb_colum, inter_group, average
            data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
            << "\t" << p_value_f << "\t" << p_value_c << "\t" << "" << "\t" << stats[2] 
            << "\t" << stats[3] << "\t" << stats[4]  << "\t" << stats[5] 
            << "\t" << stats[6] << "\t" << stats[7] << "\n";
        }

        outf.write(data.str().c_str(), data.str().size());
        itr++;
    }
}

// Quantitative Table Generation
void SnarlParser::quantitative_table(const std::vector<std::tuple<string, vector<string>, string, vector<string>>>& snarls,
                                        const std::unordered_map<std::string, double>& quantitative_phenotype, const string &chr,
                                        const std::unordered_map<std::string, std::vector<double>>& covar,
                                        const double& maf, std::vector<std::tuple<double, double, size_t>>& pvalue_vector, 
                                        const KinshipMatrix& kinship, std::ofstream& outf) {

    // Iterate over each snarl
    size_t itr = 0;
    for (const auto& tuple_snarl : snarls) {

        std::vector<std::string> list_snarl = std::get<1>(tuple_snarl);
        auto [df, allele_number] = create_quantitative_table(sampleNames, list_snarl, matrix);
        std::string snarl = std::get<0>(tuple_snarl), pos = std::get<2>(tuple_snarl);
        std::vector<std::string> type_var = std::get<3>(tuple_snarl);
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

        if (covar.size() > 0) {
            std::tuple<string, string, string, string> tuple_info = lmm_quantitative(df, quantitative_phenotype, kinship, covar);
            
            // chr, pos, snarl, type, t-dist, beta, se, p_value, allele_number
            data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
            << "\t" << std::get<0>(tuple_info) << "\t" << std::get<1>(tuple_info) 
            << "\t" << std::get<2>(tuple_info) << "\t" << std::get<3>(tuple_info) 
            << "\t" << allele_number << "\n";

        } else {

            // chr, pos, snarl, type, r2, beta, se, p_value, adjusted_pvalue, allele_number
            if (df_empty || !df_filtration) {
                pvalue_vector.push_back({1, 1, itr});
                data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
                << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA"
                << "\t" << "" << "\t" << allele_number << "\n";

            } else {
                std::tuple<string, string, string, string> tuple_info = linear_regression(df, quantitative_phenotype);
                string p_value;
                p_value = std::get<3>(tuple_info);
                pvalue_vector.push_back({string_to_pvalue(p_value), 1, itr});

                data << chr << "\t" << pos << "\t" << snarl << "\t" << type_var_str
                << "\t" << std::get<0>(tuple_info) << "\t" << std::get<1>(tuple_info) 
                << "\t" << std::get<2>(tuple_info) << "\t" << p_value << "\t" << ""
                << "\t" << allele_number << "\n";
            }
        }

        outf.write(data.str().c_str(), data.str().size());
        itr++;
    }
}

bool check_MAF_threshold_binary(const std::vector<std::vector<int>>& df, const double& maf) {

    int n = df[0].size(); // Number of columns (paths)
    int totalSum = 0;
    
    // Compute total sum of all elements in the matrix
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < n; ++j) {
            totalSum += df[i][j];
        }
    }
    
    // Check if any column's sum exceeds the threshold relative to total sum
    for (int i = 0; i < n; ++i) {
        int columnSum = df[0][i] + df[1][i];
        if (static_cast<double>(columnSum) / totalSum >= maf) {
            return false; // If any column's sum proportion meets or exceeds the threshold, return false
        }
    }

    return true; // If no column exceeds the threshold, return true
}

bool check_MAF_threshold_quantitative(const std::unordered_map<std::string, std::vector<int>>& df, const double& maf) {    
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
    const std::unordered_map<std::string, bool>& binary_groups) {

    std::unordered_map<std::string, std::vector<double>> converted_map;

    for (const auto& entry : binary_groups) {
        const std::string& sample_id = entry.first;
        bool group_value = entry.second;

        // Convert the binary value (bool) to a vector of double (1.0 or 0.0)
        std::vector<double> group_vector = {group_value ? 1.0 : 0.0};

        // Store the converted vector in the result map
        converted_map[sample_id] = group_vector;
    }

    return converted_map;
}