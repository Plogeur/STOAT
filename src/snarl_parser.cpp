#include "snarl_parser.hpp"
#include "matrix.hpp"

SnarlParser::SnarlParser(const std::string& vcf_path) : filename(vcf_path), file(vcf_path), matrix(1000000, parseHeader().size() * 2) {
    sampleNames = parseHeader();
}

std::vector<std::vector<int>> SnarlParser::create_table(
    const std::vector<std::string>& list_path_snarl)
{
    std::vector<std::vector<int>> small_grm;
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();

    // Iterate over each path_snarl in column_headers
    for (size_t idx_g = 0; idx_g < list_path_snarl.size(); ++idx_g) {
        const std::string& path_snarl = list_path_snarl[idx_g];
        const size_t number_sample = sampleNames.size();
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        std::vector<int> alleles = identify_correct_path(decomposed_snarl, row_headers_dict, matrix, number_sample*2);
        small_grm.push_back(alleles);
    }
    return small_grm;
}

std::vector<int> SnarlParser::create_table_short_path(
    const std::string& list_path_snarl)
{
    std::vector<int> small_grm;
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();

    // Iterate over each path_snarl in column_headers
    const std::string& path_snarl = list_path_snarl;
    const size_t number_sample = sampleNames.size();
    std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
    return identify_correct_path(decomposed_snarl, row_headers_dict, matrix, number_sample*2);
}

void SnarlParser::create_bim_bed(const std::unordered_map<std::string, std::vector<std::string>>& snarls,
                                  const std::string& output_bim, const std::string& output_bed) {
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
    
    int position = 0;

    // Iterate over each snarl
    for (const auto& [snarl, list_snarl] : snarls) {

        if (list_snarl.size() > 2) {continue;} // avoid multiallelic var

        // Generate a genotype table for this snarl
        std::vector<int> table = create_table_short_path(list_snarl[0]);  // 2D vector, each row = SNP, each col = sample alleles

        int chromosome = 1;  // Placeholder for chromosome number
        position += 1;  // Increment position for each SNP

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
            return std::vector<int>(num_cols, 0);
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

    // Populate allele_vector with indices of columns where all elements are 1
    std::vector<int> allele_vector;
    for (size_t col = 0; col < num_cols; ++col) {
        if (columns_all_ones[col]) {
            allele_vector.push_back(1);
        } else {
            allele_vector.push_back(0);
        }
    }
    return allele_vector;
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

void create_fam(const std::unordered_map<std::string, int>& sex, 
                const std::unordered_map<std::string, int>& pheno, 
                const std::string& output_path)

{
    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        throw std::runtime_error("Unable to open output file: " + output_path);
    }

    for (const auto& [sample, sex_code] : sex) {
        int phenotype = pheno.at(sample); // Assume pheno contains all samples.

        outfile << sample << " "       // FID (default to sample ID)
                << sample << " "       // IID (sample ID)
                << "0 0 "              // PID and MID (unknown)
                << sex_code << " "     // SEX (0 = unknown, 1 = male, 2 = female)
                << phenotype << "\n";  // PHENOTYPE (-9 = missing)
    }

    outfile.close();
}


