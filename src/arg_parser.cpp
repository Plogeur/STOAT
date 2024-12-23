#include "arg_parser.hpp"

namespace fs = std::filesystem;

// Function to parse phenotype file
void parse_pheno(const std::string &pheno_path, std::unordered_map<std::string, int> &pheno) {
    std::ifstream file(pheno_path);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << pheno_path << std::endl;
        return;
    }

    std::string line;
    bool is_header = true;

    while (std::getline(file, line)) {

        if (is_header) {
            is_header = false;
            continue;
        }
        
        std::istringstream iss(line);
        std::string fid, iid;
        int phenotype;

        // Assuming the phenotype file format: FID IID Pheno (space or tab-separated)
        if (!(iss >> fid >> iid >> phenotype)) {
            std::cerr << "Warning: Skipping malformed line: " << line << std::endl;
            continue;
        }

        // Insert into the map
        pheno[iid] = phenotype;
    }

    file.close();
}

void check_match_samples(const std::unordered_map<std::string, int>& map, const std::vector<std::string>& keys) {
    for (const auto& key : keys) {
        if (map.find(key) == map.end()) {
            throw std::runtime_error("Error: Key '" + key + "' not found in the phenotype file");
        }
    }
}

// Function to parse the snarl path file
std::unordered_map<std::string, std::vector<std::string>> parse_snarl_path(const std::string& path_file) {
    std::unordered_map<std::string, std::vector<std::string>> snarl_paths;
    std::ifstream file(path_file);

    if (!file.is_open()) {
        std::cerr << "Could not open the file: " << path_file << std::endl;
        return snarl_paths;
    }

    std::string line, snarl, path_list;
    
    // Read header
    std::getline(file, line);

    // Process each line
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::getline(ss, snarl, '\t');   // snarl column
        std::getline(ss, path_list, '\t'); // paths column
        
        if (!path_list.empty()) {
            std::istringstream path_stream(path_list);
            std::string path;
            while (std::getline(path_stream, path, ',')) {
                snarl_paths[snarl].push_back(path);
            }
        }
    }
    file.close();
    
    return snarl_paths;
}

void check_format_vcf_file(const std::string& file_path) {
    // Check if the file exists
    if (!fs::is_regular_file(file_path)) {
        throw std::invalid_argument("The file " + file_path + " does not exist.");
    }

    // Check if the file ends with .vcf
    if (file_path.size() < 4 || 
        (file_path.substr(file_path.size() - 4) != ".vcf")) {
        throw std::invalid_argument("The file " + file_path + " is not a valid VCF file. It must have a .vcf");
    }
}

std::vector<std::string> parseHeader(const std::string& file_path) {
    std::ifstream file(file_path);
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

void check_format_paths_snarl(const std::string& file_path) {
    // Check if the file exists
    if (!fs::is_regular_file(file_path)) {
        throw std::invalid_argument("The file " + file_path + " does not exist.");
    }

    // Check if the file ends with .txt or .tsv
    if (file_path.size() < 4 || 
        (file_path.substr(file_path.size() - 4) != ".txt" && file_path.substr(file_path.size() - 4) != ".tsv")) {
        throw std::invalid_argument("The file " + file_path + " is not a valid group/snarl file. It must have a .txt extension or .tsv.");
    }
}

void check_format_phenotype(const std::string& file_path) {
    // Check if the file exists
    if (!fs::is_regular_file(file_path)) {
        throw std::invalid_argument("The file " + file_path + " does not exist.");
    }

    // Open the file and check the header
    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::invalid_argument("Unable to open the file " + file_path);
    }

    std::string first_line;
    std::getline(file, first_line);
    file.close();

    std::vector<std::string> header;
    size_t start = 0;
    size_t end = first_line.find('\t');
    while (end != std::string::npos) {
        header.push_back(first_line.substr(start, end - start));
        start = end + 1;
        end = first_line.find('\t', start);
    }
    header.push_back(first_line.substr(start));

    std::vector<std::string> expected_header = {"FID", "IID", "PHENO"};
    if (header != expected_header) {
        throw std::invalid_argument("The file must contain the following headers: FID IID PHENO and be split by tabulation.");
    }
}

void parse_sex(const std::string& path, std::unordered_map<std::string, int>& sex_map) {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Unable to open sex file: " + path);
    }

    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string fid, sample_id;
        int sex_code;

        if (!(iss >> fid >> sample_id >> sex_code)) { // FID is optional for mapping
            throw std::runtime_error("Invalid format in sex file: " + line);
        }

        if (sex_code < 0 || sex_code > 2) { // Validating PLINK sex codes
            throw std::invalid_argument("Invalid sex code for sample: " + sample_id);
        }

        sex_map[sample_id] = sex_code;
    }
}
