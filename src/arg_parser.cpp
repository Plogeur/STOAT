#include "arg_parser.hpp"
#include "snarl_analyser.hpp"

namespace fs = std::filesystem;
using namespace std;

KinshipMatrix parseKinshipMatrix(const std::string& filename) {
    KinshipMatrix km;
    std::ifstream file(filename);
    std::string line;

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

const bool KinshipMatrix::empty() const {
    return ids.empty() || matrix.empty();
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

std::vector<bool> parse_binary_pheno(
    const std::string& file_path,
    const std::vector<std::string>& list_samples) {
    
    std::unordered_map<std::string, bool> binary_pheno;
    
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
            binary_pheno[iid] = static_cast<bool>(false);
        } else if (pheno == 2) {
            count_cases++;
            binary_pheno[iid] = static_cast<bool>(true);
        } else {
            throw std::runtime_error("Error: Binary phenotype must be 1 or 2");
        }
    }
    cout << "Binary phenotypes founds : " << count_controls+count_cases
    << " (Control : " << count_controls
    << ", Case : " << count_cases << ")" << endl;
    file.close();

    check_match_samples(binary_pheno, list_samples);
    std::vector<bool> vector_binary_pheno;
    vector_binary_pheno.reserve(list_samples.size());

    for (const auto& sample : list_samples) {
        auto it = binary_pheno.find(sample);
        if (it != binary_pheno.end()) {
            vector_binary_pheno.push_back(it->second);
        }
    }

    return vector_binary_pheno;
}

// Function to parse the phenotype file
std::vector<double> parse_quantitative_pheno(
    const std::string& file_path, 
    const std::vector<std::string>& list_samples) {

    std::unordered_map<std::string, double> quantitative_pheno;

    std::ifstream file(file_path);
    std::string line;
    int count_pheno = 0;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string fid, iid, phenoStr;

        if (!(iss >> fid >> iid >> phenoStr)) {
            throw std::runtime_error("Error: In parsing phenotype, malformed line: " + line);
        }

        if (firstLine) {
            firstLine = false;
            // Check that the header contains FID, IID, and PHENO
            if (fid != "FID" || iid != "IID" || phenoStr != "PHENO") {
                throw std::invalid_argument("Error: In parsing phenotype, invalid header: " + line);
            }
            continue;
        }

        try
        {
            quantitative_pheno[iid] = std::stod(phenoStr);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            throw std::runtime_error("Error: Bad phenotype type : " + phenoStr);
        }
        count_pheno++;
    }

    cout << "Quantitative phenotypes founds : " << count_pheno << endl;
    file.close();

    check_match_samples(quantitative_pheno, list_samples);
    std::vector<double> vector_quantitative_pheno;
    vector_quantitative_pheno.reserve(list_samples.size());

    for (const auto& sample : list_samples) {
        auto it = quantitative_pheno.find(sample);
        if (it != quantitative_pheno.end()) {
            vector_quantitative_pheno.push_back(it->second);
        }
    }

    return vector_quantitative_pheno;
}

// Function to open a VCF file and return pointers to the file, header, and record
std::tuple<htsFile*, bcf_hdr_t*, bcf1_t*> parse_vcf(const std::string& vcf_path) {
    // Open the VCF file
    htsFile *ptr_vcf = bcf_open(vcf_path.c_str(), "r");

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
template void check_match_samples<std::vector<double>>(const std::unordered_map<std::string, std::vector<double>>&, const std::vector<std::string>&);
template void check_match_samples<std::tuple<std::string, int, int>>(const std::unordered_map<std::string, std::tuple<std::string, int, int>>&, const std::vector<std::string>&);

template <typename T>
void check_match_samples(const std::unordered_map<std::string, T>& map, const std::vector<std::string>& keys) {
    for (const auto& key : keys) {
        if (map.find(key) == map.end()) {
            throw std::runtime_error("Error: Key '" + key + "' not found in the phenotype file");
        }
    }
    if (map.size() != keys.size()) {
        cerr << "Warning: Number of samples found in VCF does not match the number of samples in the phenotype file" << endl;
    }
}

// dict chr:string : vector{(geneName:string, sample_expression:vector<double>, start_pos:size_t, end_pos:size_t)}
std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>> 
    parse_qtl_gene_file(
    const std::string& eqtl_path, 
    const std::string& gene_position_path, 
    const std::vector<std::string>& list_samples) {

    // dict sampleName:string : vector<double> sample_expression
    auto qtl = parse_qtl_file(eqtl_path, list_samples); // and check in the same time

    // dict geneName:string : tuple{chrom:string, start_pos:size_t, end_pos:size_t}
    auto gene_position = parse_gene_positions(gene_position_path);
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::vector<double>, size_t, size_t>>> qtl_map;

    for (const auto& [gene, expression_vector] : qtl) {
        auto it = gene_position.find(gene);
        if (it != gene_position.end()) {
            const auto& [chrom, start, end] = it->second;
            qtl_map[chrom].emplace_back(gene, expression_vector, start, end);
        } else {
            std::cerr << "Error: Gene \"" << gene << "\" not found in gene positions." << std::endl;
            exit(1);
        }
    }
  
    // Warn if gene_position has more genes than qtl
    if (gene_position.size() > qtl.size()) {
        std::cerr << "Warning: More genes in the gene position file than in the QTL data." << std::endl;
    }

    return qtl_map;
}

// Function to parse the snarl path file
std::unordered_map<std::string, Snarl_data_t> parse_snarl_path(const std::string& file_path) {

    std::string line, chr, snarl, start_pos_str, end_pos_str, path_list, type_var;
    unordered_map<string, Snarl_data_t> chr_snarl_matrix;
    Snarl_data_t snarl_paths;
    std::ifstream file(file_path);
    std::string save_chr = "";

    // Read and validate header
    if (!std::getline(file, line)) {
        throw std::runtime_error("Empty file or failed to read header.");
    }
    
    // Parse actual header fields
    vector<string> header_fields;
    istringstream header_stream(line);
    string field;
    while (getline(header_stream, field, '\t')) {
        header_fields.push_back(field);
    }

    // Expected header
    vector<string> expected_header = {"CHR", "START_POS", "END_POS", "SNARL", "PATHS", "TYPE", "REF"};

    if (header_fields != expected_header) {
        // Build detailed error message
        ostringstream oss;
        oss << "Error: Invalid header format in file: " << file_path << "\n";
        oss << "  ➤ Expected: ";
        for (size_t i = 0; i < expected_header.size(); ++i) {
            oss << expected_header[i];
            if (i < expected_header.size() - 1) oss << "\\t";
        }
        oss << "\n  ➤ Got:      ";
        for (size_t i = 0; i < header_fields.size(); ++i) {
            oss << header_fields[i];
            if (i < header_fields.size() - 1) oss << "\\t";
        }
        throw runtime_error(oss.str());
    }

    // Process each line
    while (std::getline(file, line)) {
        std::istringstream ss(line);

        std::getline(ss, chr, '\t');   // chr column
        std::getline(ss, start_pos_str, '\t');   // pos column
        std::getline(ss, end_pos_str, '\t');   // pos column
        std::getline(ss, snarl, '\t');   // snarl column
        std::getline(ss, path_list, '\t'); // paths column
        std::getline(ss, type_var, '\t');   // type_var column

        std::istringstream path_stream(path_list);
        std::istringstream type_stream(type_var);
        std::vector<std::string> paths;
        std::vector<std::string> type;
        size_t start_pos = std::stoi(start_pos_str);
        size_t end_pos = std::stoi(end_pos_str);
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

        // const std::pair<size_t, size_t>& name, const Path_traversal_t& paths,
        // size_t start, size_t end, const std::vector<std::string>& path_nodes

        std::pair<size_t, size_t> snarl_pair = stringToPair(snarl);
        snarl_paths.add_snarl(snarl_pair, paths, start_pos, end_pos, type);
    }
    // last chr adding
    chr_snarl_matrix[save_chr] = std::move(snarl_paths);

    file.close();
    return chr_snarl_matrix;
}

// Function to parse the gene positions file
// dict geneName:string : tuple{chrom:string, start_pos:size_t, end_pos:size_t}
std::unordered_map<std::string, std::tuple<std::string, size_t, size_t>> parse_gene_positions(
    const std::string& filename) {

    std::unordered_map<std::string, std::tuple<std::string, size_t, size_t>> geneMap;
    std::ifstream file(filename);
    std::string line;
    
    // Read and validate header
    if (!std::getline(file, line)) {
        throw std::runtime_error("Error: Empty file or failed to read header.");
    }

    std::istringstream header_stream(line);
    std::string col1, col2, col3, col4;
    if (!(std::getline(header_stream, col1, '\t') &&
          std::getline(header_stream, col2, '\t') &&
          std::getline(header_stream, col3, '\t') &&
          std::getline(header_stream, col4, '\t')) ||
        col1 != "gene_name" || col2 != "chr" || col3 != "start" || col4 != "end") {
        throw std::runtime_error("Error: In parsing gene position file, invalid header format. Expected: gene_name\tchr\tstart\tend");
    }

    // Check for required columns
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string gene, chrom, startStr, endStr;

        std::getline(ss, gene, '\t');
        std::getline(ss, chrom, '\t');
        std::getline(ss, startStr, '\t');
        std::getline(ss, endStr, '\t');

        try {
            int start = std::stoi(startStr);
            int end = std::stoi(endStr);
            geneMap[gene] = std::make_tuple(chrom, start, end);
        } catch (...) {
            std::cerr << "Error: In parsing gene position file, invalid line: " << line << std::endl;
            exit(1);
        }
    }

    file.close();
    return geneMap;
}

// Function to parse the qtl file
// dict sampleName:string : vector<double> sample_expression
std::unordered_map<std::string, std::vector<double>> parse_qtl_file(
    const std::string& filename, const vector<std::string>& list_samples) {

    std::ifstream file(filename);
    std::unordered_map<std::string, std::vector<double>> geneExpressions;

    std::string line;
    bool isHeader = true;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        if (isHeader) {
            std::getline(ss, token, '\t');  // Skip the first column (gene name)
            std::vector<std::string> sampleNames;
            while (std::getline(ss, token, '\t')) {
                sampleNames.push_back(token);
            }

            // Check if all sample names are present in the list_samples
            for (const auto& sample : sampleNames) {
                if (std::find(list_samples.begin(), list_samples.end(), sample) == list_samples.end()) {
                    std::cerr << "Error: Sample " << sample << " not found in the list of samples." << std::endl;
                    exit(1);
                }
            }

            // warning if the number of samples in the file does not match the number of samples in the list
            if (sampleNames.size() != list_samples.size()) {
                std::cerr << "Warning: Number of samples in the qtl file is > that the number of samples in the VCF." << std::endl;
            }

            isHeader = false;  // Skip header
            continue;
        }

        std::string geneName;
        std::vector<double> expressions;

        std::getline(ss, geneName, '\t');
        while (std::getline(ss, token, '\t')) {
            try {
                expressions.push_back(std::stod(token));
            } catch (...) {
                std::cerr << "Invalid expression value for gene " << geneName << ": " << token << std::endl;
                exit(1);
            }
        }
        geneExpressions[geneName] = expressions;
    }

    file.close();
    return geneExpressions;
}

// Function to parse covariates into an unordered_map
std::vector<std::vector<double>> parse_covariates(
    const std::string& filename, 
    const std::vector<std::string>& covar_names,
    const std::vector<std::string>& list_samples) {

    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> covariate;
    std::unordered_map<string, std::vector<double>>covariate_map;

    // Read header
    std::getline(file, line);
    std::istringstream headerStream(line);
    std::vector<std::string> headers;
    std::string col;
    while (headerStream >> col) {
        headers.push_back(col);
    }

    // Check for required columns
    auto it_iid = std::find(headers.begin(), headers.end(), "IID");
    if (it_iid == headers.end()) {
        throw std::runtime_error("Error: header must include 'IID' column.\n");
        exit(1);
    }
    size_t iid_index = std::distance(headers.begin(), it_iid);

    std::unordered_map<std::string, size_t> col_index;
    for (size_t i = 0; i < headers.size(); ++i) {
        col_index[headers[i]] = i;
    }

    // Check header for covariate names
    for (const auto& name : covar_names) {
        if (col_index.find(name) == col_index.end()) {
            throw std::runtime_error("Error: covariate column '" + name + "' not found in file.\n");
            exit(1);
        }
    }

    // Read data
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::vector<std::string> tokens;
        std::string token;
        while (lineStream >> token) {
            tokens.push_back(token);
        }

        if (tokens.size() <= iid_index) continue;
        std::string iid = tokens[iid_index];

        std::vector<double> selected;
        try {
            for (const auto& name : covar_names) {
                double val = std::stod(tokens[col_index[name]]);
                selected.push_back(val);
            }
        } catch (...) {
            throw std::runtime_error("Error: Individual " + iid + " got an non-numeric value\n");
            exit(1);
        }
        covariate_map[iid] = selected;
    }

    check_match_samples(covariate_map, list_samples);

    // Order covariate_map by list_samples
    for (const auto& sample : list_samples) {
        auto it = covariate_map.find(sample);
        if (it != covariate_map.end()) {
            covariate.push_back(it->second);
        } else {
            std::cerr << "Error: Sample " << sample << " not found in the covariate file." << std::endl;
            exit(1);
        }
    }
    file.close();

    return covariate;
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
