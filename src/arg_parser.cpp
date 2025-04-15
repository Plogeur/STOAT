#include "arg_parser.hpp"
#include "snarl_parser.hpp"

namespace fs = std::filesystem;
using namespace std;

KinshipMatrix parseKinshipMatrix(const std::string& filename) {
    KinshipMatrix km;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }

    // Parse header line for IDs
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        // Skip the empty top-left cell
        std::getline(ss, token, '\t');
        while (std::getline(ss, token, '\t')) {
            km.ids.push_back(token);
        }
    }

    // Parse matrix rows
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string rowLabel;
        std::getline(ss, rowLabel, '\t'); // row label
        std::vector<double> row;
        std::string value;
        while (std::getline(ss, value, '\t')) {
            row.push_back(std::stod(value));
        }
        km.matrix.push_back(row);
    }

    file.close();
    return km;
}

std::unordered_set<std::string> parse_chromosome_reference(const string& file_path) {
    std::unordered_set<std::string> reference;
    ifstream file(file_path);
    string line;

    while (getline(file, line)) {
        reference.insert(line);
    }

    file.close();
    return reference;
}

std::unordered_map<std::string, bool> parse_binary_pheno(const std::string& file_path) {
    
    std::unordered_map<std::string, bool> parsed_pheno;
    
    std::ifstream file(file_path);
    std::string line;
    int count_controls = 0;
    int count_cases = 0;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fid, iid, phenoStr;

        if (!(iss >> fid >> iid >> phenoStr)) {
            throw std::runtime_error("Malformed line: " + line);
        }

        if (firstLine) {
            firstLine = false;
            // Check that the header contains FID, IID, and PHENO
            if (fid != "FID" || iid != "IID" || phenoStr != "PHENO") {
                throw std::invalid_argument("Invalid header: " + line);
            }
            continue;
        }

        int pheno = -1;
        try {
            pheno = std::stoi(phenoStr);
        } catch (const std::invalid_argument& e) {
            throw std::runtime_error("Bad phenotype type : " + phenoStr);
        }
        if (pheno == 1) {
            count_controls++;
            parsed_pheno[iid] = static_cast<bool>(false);
        } else if (pheno == 2) {
            count_cases++;
            parsed_pheno[iid] = static_cast<bool>(true);
        } else {
            throw std::runtime_error("Error: Phenotype must be 1 or 2");
        }
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
    int count_pheno = 0;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fid, iid, phenoStr;

        if (!(iss >> fid >> iid >> phenoStr)) {
            throw std::runtime_error("Malformed line: " + line);
        }

        if (firstLine) {
            firstLine = false;
            // Check that the header contains FID, IID, and PHENO
            if (fid != "FID" || iid != "IID" || phenoStr != "PHENO") {
                throw std::invalid_argument("Invalid header: " + line);
            }
            continue;
        }

        try
        {
            parsed_pheno[iid] = std::stod(phenoStr);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            throw std::runtime_error("Bad phenotype type : " + phenoStr);
        }
        count_pheno++;
    }

    cout << "Quantitative phenotypes founds : " << count_pheno << endl;
    file.close();
    return parsed_pheno;
}

// Function to open a VCF file and return pointers to the file, header, and record
std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path) {
    // Open the VCF file
    htsFile *ptr_vcf = bcf_open(vcf_path.c_str(), "r");
    if (!ptr_vcf) {
        throw std::runtime_error("Error: Could not open VCF file: " + vcf_path);
    }

    // Read the VCF header
    bcf_hdr_t *hdr = bcf_hdr_read(ptr_vcf);
    if (!hdr) {
        bcf_close(ptr_vcf);
        throw std::runtime_error("Error: Could not read VCF header");
    }

    // Initialize a record
    bcf1_t *rec = bcf_init();
    if (!rec) {
        bcf_hdr_destroy(hdr);
        bcf_close(ptr_vcf);
        throw std::runtime_error("Error: Failed to allocate memory for VCF record");
    }

    // Return the three initialized pointers
    return std::make_tuple(ptr_vcf, hdr, rec);
}

std::tuple<std::vector<std::string>, htsFile*, bcf_hdr_t*, bcf1_t*> parseHeader(const std::string& vcf_path) {
    auto [ptr_vcf, hdr, rec] = parse_vcf(vcf_path);

    std::vector<std::string> list_samples;
    // Get the samples names
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        list_samples.push_back(bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i));
    }
        
    return std::make_tuple(list_samples, ptr_vcf, hdr, rec);
}

// Explicit instantiation for specific types
template void check_match_samples<bool>(const std::unordered_map<std::string, bool>&, const std::vector<std::string>&);
template void check_match_samples<double>(const std::unordered_map<std::string, double>&, const std::vector<std::string>&);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys) {
    for (const auto& key : keys) {
        if (map.find(key) == map.end()) {
            throw std::runtime_error("Error: Key '" + key + "' not found in the phenotype file");
        }
    }
    if (map.size() != keys.size()) {
        cerr << "Warning:  Number of samples found in VCF does not match the number of samples in the phenotype file" << endl;
    }
}

// Function to parse the snarl path file
std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> parse_snarl_path(const std::string& file_path) {

    std::string line, chr, pos, snarl, path_list, type_var;
    unordered_map<string, std::vector<std::tuple<string, vector<string>, string, vector<string>>>> chr_snarl_matrix;
    std::vector<std::tuple<string, vector<string>, string, vector<string>>> snarl_paths;
    std::ifstream file(file_path);
    std::string save_chr = "";

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
        int size_paths = 0;

        // create a vector of paths
        while (std::getline(path_stream, path_list, ',')) {
            size_paths++;
            paths.push_back(path_list);
        }

        // create a vector of types
        while (std::getline(type_stream, type_var, ',')) {
            type.push_back(type_var);
        }

        if (chr != save_chr && !save_chr.empty()) {
            chr_snarl_matrix[save_chr] = std::move(snarl_paths);
            snarl_paths.clear();
        }
        save_chr = chr;

        // {snarl, paths, chr, pos, type} 
        snarl_paths.push_back(make_tuple(snarl, paths, pos, type));
    }
    // last chr adding
    chr_snarl_matrix[save_chr] = std::move(snarl_paths);

    file.close();
    return chr_snarl_matrix;
}

QTL parseExpressionFile(const std::string& filename) {
    QTL data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // Parse header line for sample IDs
    if (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        // First column is "Gene_ID", skip it
        std::getline(ss, token, '\t');

        // Sample IDs
        while (std::getline(ss, token, '\t')) {
            data.sample_ids.push_back(token);
        }
    }

    // Parse gene expression rows
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string gene_id;
        std::string value;
        std::vector<double> expression_values;

        // First column: gene ID
        std::getline(ss, gene_id, '\t');
        data.gene_ids.push_back(gene_id);

        // Remaining columns: expression values
        while (std::getline(ss, value, '\t')) {
            expression_values.push_back(std::stod(value));
        }

        data.expression_matrix.push_back(expression_values);
    }

    file.close();
    return data;
}

// Function to check covariate format
void check_format_covariate(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
    }

    std::string line;
    int numCols = -1; // Store column count
    int lineCount = 0;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        int colCount = 0;

        while (ss >> token) colCount++; // Count columns in this line

        if (numCols == -1) {
            numCols = colCount; // Set number of columns from first line
        } else if (colCount != numCols) {
            std::cerr << "Error: Inconsistent column count in row " << lineCount + 1 << std::endl;
        }
        lineCount++;
    }
    file.close();
}

// Function to parse covariates into an unordered_map
std::unordered_map<std::string, std::vector<double>> parseCovariate(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return {};
    }

    std::unordered_map<std::string, std::vector<double>> covariates;
    std::string line;
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string sample_id;
        ss >> sample_id;

        std::vector<double> values;
        double value;
        while (ss >> value) {
            values.push_back(value);
        }

        covariates[sample_id] = values;
    }
    
    file.close();
    return covariates;
}

template void check_phenotype_covariate<double>(const std::unordered_map<std::string, double>& phenotype, 
    const std::unordered_map<std::string, std::vector<double>>& covariates);

template void check_phenotype_covariate<bool>(const std::unordered_map<std::string, bool>& phenotype, 
    const std::unordered_map<std::string, std::vector<double>>& covariates);

// Function to check if phenotype and covariate files match
template <typename T>
void check_phenotype_covariate(const std::unordered_map<std::string, T>& phenotype, 
    const std::unordered_map<std::string, std::vector<double>>& covariates) {

        // Check if all phenotype samples are present in covariates
        for (const auto& pair : phenotype) {
        if (covariates.find(pair.first) == covariates.end()) {
            std::cerr << "Error: Missing covariate data for sample " << pair.first << std::endl;
            EXIT_FAILURE;
        }
    }
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
