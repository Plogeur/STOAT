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

std::string set_precision_chi2(const boost::multiprecision::cpp_dec_float_50& value) {
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

// Zaykin et al. generic p-value combination
double combine_pvalue_from_strings(const std::string& p1, const std::string& p2) {
    bool na1 = is_na(p1);
    bool na2 = is_na(p2);

    if (!na1 && !na2) {
        double d1 = std::stod(p1);
        double d2 = std::stod(p2);

        // Apply Fisher's transformation
        double T1 = -2.0 * std::log(d1);
        double T2 = -2.0 * std::log(d2);
        double Y = T1 + T2;

        // Degrees of freedom = 4 (2 tests)
        boost::math::chi_squared dist(4);
        return 1.0 - boost::math::cdf(dist, Y);

    } else if (!na1) {
        return std::stod(p1);
    } else if (!na2) {
        return std::stod(p2);
    } else {
        return 1.0;  // both are NA
    }
}

double string_to_pvalue(const std::string& p1) {
    bool na1 = is_na(p1);

    if (!na1) {
        return std::stod(p1);
    } else {
        return 1.0;
    }
}

// Function to check significance from a string
bool isPValueSignificant(size_t numDigits, const std::string& pvalue_str) {
    double pvalue;
    try {
        if (pvalue_str == "NA") {
            return false; // Treat "NA" as not significant
        } else {
            pvalue = std::stod(pvalue_str);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error parsing pvalue string : " << pvalue_str << " " << e.what() << "\n";
        return false;
    }

    double threshold = std::pow(10.0, -static_cast<int>(numDigits));
    return pvalue < threshold;
}

// Write the table to a TSV file
void writeSignificantTableToTSV(
    const std::vector<std::vector<size_t>>& table,
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
