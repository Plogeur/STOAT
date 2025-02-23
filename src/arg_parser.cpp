#include "arg_parser.hpp"

namespace fs = std::filesystem;
using namespace std;

std::unordered_map<std::string, bool> parse_binary_pheno(const std::string& file_path) {
    
    std::unordered_map<std::string, bool> parsed_pheno;
    std::ifstream file(file_path);

    std::string line;
    std::getline(file, line);
    int count_controls = 0;
    int count_cases = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fid, iid, phenoStr;

        if (!(iss >> fid >> iid >> phenoStr)) {
            throw std::runtime_error("Malformed line: " + line);
        }
        int pheno = std::stoi(phenoStr);
        if (pheno == 0) {
            count_controls++;
        } else if (pheno == 1) {
            count_cases++;
        } else {
            throw std::runtime_error("Error: Phenotype must be 0 or 1");
        }
        parsed_pheno[iid] = static_cast<bool>(pheno);
    }
    cout << "Binary phenotypes founds : " << count_controls+count_cases
    << " (Control : " << count_controls
    << ", Case : " << count_cases << ")" << endl;
    file.close();
    return parsed_pheno;
}

// Function to parse the phenotype file
std::unordered_map<std::string, double> parse_quantitative_pheno(const std::string& file_path) {
    std::unordered_map<std::string, double> parsed_pheno;

    std::ifstream file(file_path);
    std::string line;
    std::getline(file, line);  // Skip header line
    int count_pheno = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string col1, iid;
        double pheno;

        // Read the first column, IID (second column), and PHENO (third column)
        std::getline(iss, col1, '\t');  // Skip the first column
        std::getline(iss, iid, '\t');   // Get the IID
        iss >> pheno;                   // Read the PHENO as a double

        parsed_pheno[iid] = pheno;  // Add to map
        count_pheno++;
    }

    cout << "Quantitative phenotypes founds : " << count_pheno << endl;
    file.close();
    return parsed_pheno;
}

std::vector<std::string> parseHeader(const std::string& vcf_path) {

    // Open the VCF file
    htsFile *vcf_file = bcf_open(vcf_path.c_str(), "r");
    if (!vcf_file) {
        throw std::runtime_error("Error: Could not open VCF file: " + vcf_path);
    }
    
    // Read the VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(vcf_file);
    if (!hdr) {
        bcf_close(vcf_file);
        throw std::runtime_error("Error: Could not read VCF header");
    }

    std::vector<std::string> sampleNames;

    // Get the samples names
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        sampleNames.push_back(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }

    // Cleanup
    bcf_hdr_destroy(hdr);
    bcf_close(vcf_file);
    return sampleNames;
}

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys) {
    for (const auto& key : keys) {
        if (map.find(key) == map.end()) {
            throw std::runtime_error("Error: Key '" + key + "' not found in the phenotype file");
        }
    }
    if (map.size() != keys.size()) {
        cerr << "Warning:  Number of samples found in VCF does not match the number of samples in the phenotype file" << endl;
        cerr << "Number of samples found in VCF: " << keys.size() << endl;
        cerr << "Number of samples in the phenotype file: " << map.size() << endl;
    }
}

// Explicit instantiation for specific types
template void check_match_samples<bool>(const std::unordered_map<std::string, bool>&, const std::vector<std::string>&);
template void check_match_samples<double>(const std::unordered_map<std::string, double>&, const std::vector<std::string>&);

// Function to parse the snarl path file
std::vector<std::tuple<string, vector<string>, string, string, vector<string>>> parse_snarl_path(const std::string& file_path) {

    std::string line, chr, pos, snarl, path_list, type_var;
    std::vector<std::tuple<string, vector<string>, string, string, vector<string>>> snarl_paths;
    std::ifstream file(file_path);

    // Read header
    std::getline(file, line);

    // Process each line
    while (std::getline(file, line)) {
        std::istringstream ss(line);

        std::getline(ss, chr, '\t');   // chr column
        std::getline(ss, pos, '\t');   // pos column
        std::getline(ss, snarl, '\t');   // snarl column
        std::getline(ss, path_list, '\t'); // paths column
        std::getline(ss, type_var, '\t');   // type_var column

        std::istringstream path_stream(path_list);
        std::istringstream type_stream(type_var);
        std::vector<std::string> paths;
        std::vector<std::string> type;

        // create a vector of paths
        while (std::getline(path_stream, path_list, ',')) {
            paths.push_back(path_list);
        }

        // create a vector of types
        while (std::getline(type_stream, type_var, ',')) {
            type.push_back(type_var);
        }

        // {snarl, paths, chr, pos, type} 
        snarl_paths.push_back(make_tuple(snarl, paths, chr, pos, type));

    }

    file.close();
    return snarl_paths;
}

void check_format_quantitative_phenotype(const std::string& file_path) {
    
    std::ifstream file(file_path);
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
        throw std::invalid_argument("The file must contain the following headers: FID, IID, PHENO and be split by tabulation.");
    }
}

bool is_valid_iid(const std::string& iid) {
    // Check if IID is valid (e.g., not empty, no special characters)
    return !iid.empty();
}

bool is_valid_group(bool group) {
    // Check if group is either 0 or 1 (valid boolean)
    return group == 0 || group == 1;
}

void check_file(const std::string& file_path) {
    
    if (!fs::is_regular_file(file_path)) {
        throw std::invalid_argument("The file " + file_path + " does not exist.");
    }

    std::ifstream file(file_path);
    if (!file.is_open()) {
        throw std::invalid_argument("Unable to open the file " + file_path);
    }

    file.close();
}

void check_format_binary_phenotype(const std::string& file_path) {

    std::ifstream file(file_path);

    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Error opening file: " + file_path);
    }

    std::string line;
    bool firstLine = true;

    while (std::getline(file, line)) {

        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string fid, iid, pheno;

        // Parse FID, IID, and PHENO from the line
        if (!std::getline(ss, fid, '\t') || 
            !std::getline(ss, iid, '\t') || 
            !std::getline(ss, pheno, '\t')) {
            throw std::invalid_argument("Invalid line format: " + line);
        }

        if (firstLine) {
            firstLine = false;
            // Check that the header contains FID, IID, and PHENO
            if (fid != "FID" || iid != "IID" || pheno != "PHENO") {
                throw std::invalid_argument("Invalid header: " + line);
            }
            continue;
        }

        // Check if PHENO is 0 or 1
        if (pheno != "0" && pheno != "1") {
            throw std::invalid_argument("Invalid PHENO value: " + pheno + " in line: " + line);
        }
    }

    // If no errors occurred, the file is valid
    file.close();
}
