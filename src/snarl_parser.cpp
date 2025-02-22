#include "snarl_parser.hpp"
#include "matrix.hpp"
#include "binary_analysis.hpp"
#include "quantitative_analysis.hpp"

using namespace std;

SnarlParser::SnarlParser(const std::string& vcf_path) : filename(vcf_path), file(vcf_path) {
    sampleNames = parseHeader();
    matrix = Matrix(1000000, sampleNames.size() * 2);
}

std::vector<std::string> SnarlParser::parseHeader() {

    std::vector<std::string> sampleNames;
    
    return sampleNames;
}

// Function to extract an integer from a string starting at index `i`
std::pair<int, std::string> determine_str(const std::string& s, int length_s, int i) {
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
void SnarlParser::pushMatrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn) {
    // Retrieve or add the index in one step and calculate length once
    size_t lengthOrderedMap = rowHeaderDict.size();
    size_t idxSnarl = getOrAddIndex(rowHeaderDict, decomposedSnarl, lengthOrderedMap);

    // Check if a new matrix chunk is needed
    size_t currentRowsNumber = matrix.getRows();
    if (lengthOrderedMap > currentRowsNumber - 1) {
        matrix.expandMatrix();
    }

    // Add data to the matrix
    matrix.set(idxSnarl, indexColumn);
}
 // Function to parse VCF and process genotypes
void parse_vcf(const std::string &vcf_path, Matrix &matrix) {
    htsFile *vcf_file = bcf_open(vcf_path.c_str(), "r");
    if (!vcf_file) {
        std::cerr << "Error: Could not open VCF file " << vcf_path << "\n";
        return;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcf_file);
    if (!hdr) {
        std::cerr << "Error: Could not read VCF header\n";
        bcf_close(vcf_file);
        return;
    }

    bcf1_t *rec = bcf_init();
    if (!rec) {
        std::cerr << "Error: Failed to allocate memory for record\n";
        bcf_hdr_destroy(hdr);
        bcf_close(vcf_file);
        return;
    }

    std::map<int, std::string> row_header_dict;

    while (bcf_read(vcf_file, hdr, rec) >= 0) {
        bcf_unpack(rec, BCF_UN_STR);

        // Skip variant if LV != 0
        int32_t *lv_value = NULL;
        int n_lv = 0;
        n_lv = bcf_get_info_int32(hdr, rec, "LV", &lv_value, &n_lv);
        if (n_lv > 0 && lv_value[0] != 0) {
            free(lv_value);
            continue;
        }
        free(lv_value);

        // Extract AT field and split by ','
        int32_t *at_values = NULL;
        int nat = 0;
        nat = bcf_get_info_int32(hdr, rec, "AT", &at_values, &nat);
        std::vector<std::string> snarl_list;
        if (nat > 0 && at_values) {
            for (int i = 0; i < nat; i++) {
                snarl_list.push_back(std::to_string(at_values[i]));
            }
        }
        free(at_values);

        // Decompose snarl list
        std::vector<std::vector<sstring>> list_list_decomposed_snarl = decompose_snarl(snarl_list);

        // Extract GT field
        int ngt = 0;
        int32_t *gt = NULL;
        ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt);

        if (ngt > 0 && gt != NULL) {
            for (int index_column = 0; index_column < rec->n_sample; index_column++) {
                int allele_1 = bcf_gt_allele(gt[index_column * 2]);
                int allele_2 = bcf_gt_allele(gt[index_column * 2 + 1]);
                int col_idx = index_column * 2;

                if (allele_1 != -1) {
                    for (int decompose_allele_1 : list_list_decomposed_snarl[allele_1]) {
                        matrix.push_matrix(allele_1, decompose_allele_1, row_header_dict, col_idx);
                    }
                }

                if (allele_2 != -1) {
                    for (int decompose_allele_2 : list_list_decomposed_snarl[allele_2]) {
                        matrix.push_matrix(allele_2, decompose_allele_2, row_header_dict, col_idx + 1);
                    }
                }
            }
        }
        free(gt);
    }

    // Set row headers
    matrix.set_row_header(row_header_dict);

    // Cleanup
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf_file);
}

std::vector<int> identify_correct_path(
    const std::vector<std::string>& decomposed_snarl,
    const std::unordered_map<std::string, size_t>& row_headers_dict,
    const Matrix& matrix,
    const size_t num_cols
) {
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
void SnarlParser::binary_table(const unordered_map<string, tuple<vector<string>, string, string, vector<string>>>& snarls,
                                  const std::unordered_map<std::string, bool>& binary_groups,
                                  const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);

    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tP_FISHER\tP_CHI2\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, tuple_snarl] : snarls) {

        std::vector<std::string> list_snarl = std::get<0>(tuple_snarl);
        std::vector<std::vector<int>> df = create_binary_table(binary_groups, list_snarl, sampleNames, matrix);
        std::vector<std::string> stats = binary_stat_test(df);

        std::string chrom = std::get<1>(tuple_snarl), pos = std::get<2>(tuple_snarl);
        std::vector<std::string> type_var = std::get<3>(tuple_snarl);

        // make a string separated by , from a vector of string
        std::ostringstream oss;
            for (size_t i = 0; i < type_var.size(); ++i) {
                if (i != 0) oss << ","; // Add comma before all elements except the first
                oss << type_var[i];
            }

        std::string type_var_str = oss.str();
        
        // fisher_p_value, chi2_p_value 
        // TODO : add other metrics 
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var_str << "\t"
             << "\t" << stats[0] << "\t" << stats[1] << "\n";
        
        outf.write(data.str().c_str(), data.str().size());
    }
}

// Quantitative Table Generation
void SnarlParser::quantitative_table(const unordered_map<string, tuple<vector<string>, string, string, vector<string>>>& snarls,
                                        const std::unordered_map<std::string, double>& quantitative_phenotype,
                                        const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);
    
    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tSE\tBETA\tP\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, tuple_snarl] : snarls) {

        std::vector<std::string> list_snarl = std::get<0>(tuple_snarl);
        std::unordered_map<std::string, std::vector<int>> df = create_quantitative_table(sampleNames, list_snarl, matrix);

        // std::make_tuple(se, beta, p_value)
        std::tuple<std::string, std::string, std::string> tuple_info = linear_regression(df, quantitative_phenotype);

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt
            << "\t" << std::get<0>(tuple_info) << "\t" << std::get<1>(tuple_info) 
            << "\t" << std::get<2>(tuple_info) << "\n";

        outf.write(data.str().c_str(), data.str().size());
    }
}

// Extracts the genotype from the genotype string
std::vector<int> extractGenotype(const std::string& genotypeStr) {
    std::vector<int> alleles;
    std::vector<std::string> alleleStrs = split(genotypeStr, '/');

    for (const std::string& alleleStr : alleleStrs) {
        if (alleleStr[0] == '.') {
            alleles.push_back(-1); // Use -1 for missing data
        } else {
            alleles.push_back(std::stoi(alleleStr)); // Convert allele string to integer
        }
    }
    return alleles;
}

// Extracts the AT field from the INFO field and returns a vector of strings
std::vector<std::string> extractATField(const std::string& infoField) {
    std::vector<std::string> atValues; // Change to vector of strings
    std::vector<std::string> infoParts = split(infoField, ';');
    
    for (const std::string& part : infoParts) {
        if (part.rfind("AT=", 0) == 0) {
            std::string atData = part.substr(3);  // Remove "AT="
            // Split the AT data using ',' only
            atValues = split(atData, ','); // Directly assign the split result to atValues
            return atValues; // Return the result
        }
    }
    return atValues;  // Return empty if AT field is not found
}

// Utility function to split strings
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}
