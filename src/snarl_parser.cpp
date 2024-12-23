#include "snarl_parser.hpp"
#include "matrix.hpp"
#include "binary_analysis.hpp"
#include "quantitative_analysis.hpp"

SnarlParser::SnarlParser(const std::string& vcf_path) : filename(vcf_path), file(vcf_path), matrix(1000000, parseHeader().size() * 2) {
    sampleNames = parseHeader();
}

std::vector<std::string> SnarlParser::parseHeader() {
    file.clear();                     // Clear any flags in the file stream
    file.seekg(0, std::ios::beg);     // Move the file stream to the beginning

    std::vector<std::string> sampleNames;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;                 // Skip any empty lines
        }

        // Stop reading the header once we encounter a non-comment line
        if (line[0] != '#') {
            file.seekg(-static_cast<int>(line.size()) - 1, std::ios::cur);  // Go back to the start of this line
            break;
        }

        // Process the header line that starts with "#CHROM"
        if (line.substr(0, 6) == "#CHROM") {
            std::istringstream headerStream(line);
            std::string sampleName;

            // Skip the mandatory VCF columns
            for (int i = 0; i < 9; ++i) {
                headerStream >> sampleName;
            }

            // Read remaining entries as sample names
            while (headerStream >> sampleName) {
                sampleNames.push_back(sampleName);
            }
        }
    }
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

// Main function that parses the VCF file and fills the matrix
void SnarlParser::fill_matrix() {

    std::unordered_map<std::string, size_t> row_header_dict;

    for (std::string line; std::getline(file, line);) {

        Variant variant(line, sampleNames);
        const std::vector<std::vector<int>>& genotypes = variant.genotypes;
        const std::vector<std::string>& snarl_list = variant.atInfo;
        const std::vector<std::vector<std::string>> list_list_decomposed_snarl = decompose_snarl(snarl_list);

        // Process genotypes and fill matrix
        for (size_t index_column = 0; index_column < genotypes.size(); ++index_column) {

            int allele_1 = genotypes[index_column][0];
            int allele_2 = genotypes[index_column][1];
            size_t col_idx = index_column * 2;

            if (allele_1 == -1 || allele_2 == -1) { // Handle missing genotypes (./.)
                continue;
            }

            // Measure the time taken to push each allele into the matrix
            for (auto& decompose_allele_1 : list_list_decomposed_snarl[allele_1]) {
                pushMatrix(decompose_allele_1, row_header_dict, col_idx);
            }

            for (auto& decompose_allele_2 : list_list_decomposed_snarl[allele_2]) {
                pushMatrix(decompose_allele_2, row_header_dict, col_idx+1);
            }
        }
        matrix.set_row_header(row_header_dict);
    }
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
void SnarlParser::binary_table(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
                                  const std::unordered_map<std::string, bool>& binary_groups,
                                  const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);

    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tP_Fisher\tP_Chi2\tGROUP_1_PATH_1\tGROUP_1_PATH_2\tGROUP_2_PATH_1\tGROUP_2_PATH_2\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, list_snarl] : snarls) {

        std::vector<std::vector<int>> df = create_binary_table(binary_groups, list_snarl, sampleNames, matrix);
        std::vector<std::string> stats = binary_stat_test(df);

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        // Stats is a vector containing the values in the order: 
        // fisher_p_value, chi2_p_value, GIPI, GIPII, GIIPI, GIIPII
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt
             << "\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2] << "\t" << stats[3]
             << "\t" << stats[4] << "\t" << stats[5] << "\n";
        
        outf.write(data.str().c_str(), data.str().size());
    }
}

// Quantitative Table Generation
void SnarlParser::quantitative_table(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
                                        const std::unordered_map<std::string, float>& quantitative_phenotype,
                                        const std::string& output) 
{
    std::ofstream outf(output, std::ios::binary);
    
    // Write headers
    std::string headers = "CHR\tPOS\tSNARL\tTYPE\tREF\tALT\tSE\tBETA\tP\n";
    outf.write(headers.c_str(), headers.size());

    // Iterate over each snarl
    for (const auto& [snarl, list_snarl] : snarls) {

        std::unordered_map<std::string, std::vector<int>> df = create_quantitative_table(list_snarl, sampleNames, matrix);

        // std::make_tuple(se, beta, p_value)
        std::tuple<double, double, double> tuple_info = linear_regression(df, quantitative_phenotype);

        std::string chrom = "NA", pos = "NA", type_var = "NA", ref = "NA", alt = "NA";
        std::stringstream data;
        data << chrom << "\t" << pos << "\t" << snarl << "\t" << type_var << "\t" << ref << "\t" << alt
            << "\t" << std::get<0>(tuple_info) << "\t" << std::get<1>(tuple_info) 
            << "\t" << std::get<2>(tuple_info) << "\n";

        outf.write(data.str().c_str(), data.str().size());
    }
}

// Parses a variant line from the VCF file and extracts genotype and AT field
Variant::Variant(std::string& line, std::vector<std::string> sampleNames) {
    std::istringstream variantStream(line);
    std::string columnData;
    std::string infoField;
    std::string genotypeStr;

    // Skip the first 9 columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
    for (int i = 0; i < 9; ++i) {
        variantStream >> columnData;
        if (i == 7) {
            infoField = columnData; // Capture the INFO column for parsing
        }
    }

    // Read the genotype data for each sample
    for (const auto& sample : sampleNames) {
        if (!(variantStream >> genotypeStr)) {
            throw std::runtime_error("Error reading genotype for sample: " + sample);
        }
        genotypes.push_back(extractGenotype(genotypeStr));
    }
    // Extract AT field from the INFO column
    atInfo = extractATField(infoField);
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
