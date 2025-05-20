#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

// Type aliases
using FeatureMatrix = std::vector<std::vector<double>>;
using PhenotypeVector = std::vector<double>;

// Function to parse the feature file
void parse_feature_file(
    const std::string& feature_filename,
    std::vector<std::string>& sample_ids,
    FeatureMatrix& features
) {
    std::ifstream infile(feature_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open feature file");
    }

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;

        std::string sample_id;
        std::getline(ss, sample_id, '\t');
        sample_ids.push_back(sample_id);

        std::vector<double> feature_row;
        while (std::getline(ss, token, '\t')) {
            feature_row.push_back(std::stod(token));
        }

        features.push_back(feature_row);
    }
}

// Function to parse the phenotype file
void parse_phenotype_file(
    const std::string& phenotype_filename,
    const std::vector<std::string>& sample_ids,
    PhenotypeVector& phenotype
) {
    std::ifstream infile(phenotype_filename);
    if (!infile) {
        throw std::runtime_error("Unable to open phenotype file");
    }

    std::unordered_map<std::string, double> phenotype_map;

    std::string line;
    std::getline(infile, line);  // skip header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string fid, iid, pheno_str;
        std::getline(ss, fid, '\t');
        std::getline(ss, iid, '\t');
        std::getline(ss, pheno_str, '\t');

        phenotype_map[iid] = std::stod(pheno_str);
    }

    for (const auto& sample : sample_ids) {
        if (phenotype_map.find(sample) != phenotype_map.end()) {
            phenotype.push_back(phenotype_map[sample]);
        } else {
            throw std::runtime_error("Sample ID not found in phenotype file: " + sample);
        }
    }
}

// Example usage
int main() {
    std::string feature_file = "../data/quantitative/phenotype.tsv";
    std::string phenotype_file = "../output/regression/4220_4223.tsv";

    std::vector<std::string> sample_ids;
    FeatureMatrix features;
    PhenotypeVector phenotype;

    try {
        parse_feature_file(feature_file, sample_ids, features);
        parse_phenotype_file(phenotype_file, sample_ids, phenotype);

        std::cout << "Parsed " << features.size() << " samples with "
                  << features[0].size() << " features.\n";
        std::cout << "Parsed " << phenotype.size() << " phenotype values.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return 0;
}


// g++ -std=c++17 -o parsing_regression parsing_regression.cpp
// g++ -std=c++17 -o parsing_regression parsing_regression.cpp