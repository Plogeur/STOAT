#include "adjusted_pvalue.hpp"
#include "utils.hpp"

// pvalue, adjusted_pvalue, snarl_index
void adjust_pvalues_BH(
    std::vector<std::tuple<double, double, size_t>>& data) {
    size_t n = data.size();

    // Sort by raw p-value (ascending)
    std::sort(data.begin(), data.end(),
              [](const auto& a, const auto& b) {
                  return std::get<0>(a) < std::get<0>(b);
              });

    // Apply Benjamini-Hochberg procedure
    std::vector<double> adjusted(n);
    for (size_t i = 0; i < n; ++i) {
        double p = std::get<0>(data[i]);
        adjusted[i] = p * n / (i + 1);  // BH formula
    }

    // Ensure monotonicity (non-decreasing order after adjustment)
    for (size_t i = n - 1; i > 0; --i) {
        if (adjusted[i - 1] > adjusted[i]) {
            adjusted[i - 1] = adjusted[i];
        }
    }

    // Clamp to [0,1] and update adjusted p-values in tuples
    for (size_t i = 0; i < n; ++i) {
        std::get<1>(data[i]) = std::min(1.0, adjusted[i]);
    }

    // Sort back to original index order
    std::sort(data.begin(), data.end(),
              [](const auto& a, const auto& b) {
                  return std::get<2>(a) < std::get<2>(b);
              });
}

void update_gwas_file_with_adjusted_pvalues(
    const std::string& input_file,
    const std::string& output_file_significant,
    const int column_index,
    const std::vector<std::tuple<double, double, size_t>>& adjusted_pvalues) {

    std::ifstream infile(input_file);
    std::ofstream temp_file("temp_gwas.tsv");  // Temporary file to overwrite original
    std::ofstream outfile_significant(output_file_significant);

    std::string line;
    size_t line_index = 0;

    // Process header
    if (std::getline(infile, line)) {
        temp_file << line << '\n';
        outfile_significant << line << '\n';
    }

    // Process each data line
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        // Update the adjusted p-value column
        double adjusted_p = std::get<1>(adjusted_pvalues[line_index]);
        columns[column_index] = set_precision(adjusted_p);

        // Reconstruct the line
        std::ostringstream new_line;
        for (size_t i = 0; i < columns.size(); ++i) {
            new_line << columns[i];
            if (i != columns.size() - 1) new_line << '\t';
        }

        // Write updated line to temp file
        temp_file << new_line.str() << '\n';

        // If significant, write to separate file
        if (adjusted_p < 1e-5) {
            outfile_significant << new_line.str() << '\n';
        }

        ++line_index;
    }

    infile.close();
    temp_file.close();
    outfile_significant.close();

    // Replace the original file with the updated file
    std::remove(input_file.c_str());
    std::rename("temp_gwas.tsv", input_file.c_str());
}
