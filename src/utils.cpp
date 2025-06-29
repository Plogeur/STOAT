#include "utils.hpp"

std::string set_precision(const double& value) {
    std::ostringstream oss;
    oss << std::setprecision(4);

    if (value < 1e-4) {
        oss << std::scientific << value; // Scientific notation with 4 decimals
    } else {
        oss << std::fixed << value; // Fixed-point notation with 4 decimals
    }
    return oss.str();
}

std::string set_precision_float_50(const boost::multiprecision::cpp_dec_float_50& value) {
    std::ostringstream oss;
    oss << std::setprecision(4);
    if (value < boost::multiprecision::cpp_dec_float_50("1e-4")) {
        oss << std::scientific << value;
    } else {
        oss << std::fixed << value;
    }
    return oss.str();
}

bool is_na(const std::string& s) {
    return s.empty() || s == "NA";
}

double string_to_pvalue(const std::string& p1) {
    bool na1 = is_na(p1);

    if (!na1) {
        return std::stod(p1);
    } else {
        return 1.0;
    }
}

// Function to check significance from a std::string
bool isPValueSignificant(const double& pvalue_threshold, const std::string& pvalue_str) {
    double pvalue;
    try {
        if (pvalue_str == "NA") {
            return false; // Treat "NA" as not significant
        } else {
            pvalue = std::stod(pvalue_str);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing pvalue std::string : " << pvalue_str << " " << e.what() << "\n";
        return false;
    }
    return pvalue < pvalue_threshold;
}

// Write the table to a TSV file
void writeSignificantTableToTSV(
    const std::vector<std::vector<double>>& table,
    const std::vector<std::string>& list_snarl,
    const std::vector<std::string>& list_samples,
    const std::string& filename) {

    std::ofstream outFile(filename);

    // Write header
    outFile << "sample_name";
    for (const auto& snarl_name : list_snarl) {
        outFile << "\t" << snarl_name;
    }
    outFile << "\n";

    // Write each sample's data
    size_t itr = 0;
    for (const auto& allele_vector : table) {
        outFile << list_samples[itr];

        for (size_t i=0; i < allele_vector.size(); ++i) {
            outFile << "\t" << allele_vector[i];
        }
        outFile << "\n";
        ++itr;
    }
    outFile.close();
}

// Adjust p-values using Holm-Bonferroni correction
std::vector<double> adjusted_holm(const std::vector<double>& p_values) {
    int m = p_values.size();
    std::vector<std::pair<double, int>> indexed;
    for (int i = 0; i < m; ++i) {
        indexed.emplace_back(p_values[i], i);
    }

    // Sort by p-value
    std::sort(indexed.begin(), indexed.end());

    std::vector<double> adjusted(m);
    double prev = 0.0;
    for (int i = 0; i < m; ++i) {
        double raw = (m - i) * indexed[i].first;
        raw = std::min(raw, 1.0);
        adjusted[i] = std::max(prev, raw); // ensure monotonicity
        prev = adjusted[i];
    }

    // Reorder to original positions
    std::vector<double> reordered(m);
    for (int i = 0; i < m; ++i) {
        reordered[indexed[i].second] = adjusted[i];
    }

    return reordered;
}

void retain_indices(std::vector<double>& vec, const std::unordered_set<size_t>& indices_to_keep) {
    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx < vec.size(); ++read_idx) {
        if (indices_to_keep.count(read_idx)) {
            vec[write_idx++] = vec[read_idx];
        }
    }
    vec.resize(write_idx);
}

template std::string vectorToString(const std::vector<std::string>& vec);
template std::string vectorToString(const std::vector<size_t>& vec);

template<typename T>
std::string vectorToString(const std::vector<T>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec[i];
    }
    return oss.str();
}

template std::vector<std::string> stringToVector(const std::string& vec);
template std::vector<size_t> stringToVector(const std::string& vec);

template <typename T>
std::vector<T> stringToVector(const std::string& str) {
    std::vector<T> result;
    std::istringstream iss(str);
    std::string token;

    while (std::getline(iss, token, ',')) {
        std::istringstream tokenStream(token);
        T value;
        tokenStream >> value;
        if (tokenStream.fail()) {
            throw std::runtime_error("Failed to parse token: " + token);
        }
        result.push_back(value);
    }

    return result;
}

std::string vectorPathToString(const std::vector<Path_traversal_t>& vec) {
    std::ostringstream oss;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) oss << ",";
        oss << vec[i].to_string();
    }
    return oss.str();
}

std::string get_sample_name_from_path(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path) {

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        return graph.get_locus_name(path);
    } else {
        return graph.get_sample_name(path);
    }

}

sample_hap_t get_sample_and_haplotype(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path) {
    sample_hap_t result;

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        result.sample = graph.get_locus_name(path);
    } else {
        result.sample = graph.get_sample_name(path);
    }
    result.haplotype = graph.get_haplotype(path);

    return result;
}
