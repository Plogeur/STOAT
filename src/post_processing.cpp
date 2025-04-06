#include "post_processing.hpp"
#include "utils.hpp"

// BH Adjustment
void adjust_pvalues_BH(std::vector<std::tuple<double, double, size_t>>& data) {
    size_t n = data.size();
    std::sort(data.begin(), data.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

    std::vector<double> adjusted(n);
    for (size_t i = 0; i < n; ++i) {
        double p = std::get<0>(data[i]);
        if (p == 0.0) continue;
        adjusted[i] = p * n / (i + 1);
    }

    for (size_t i = n - 1; i > 0; --i) {
        if (adjusted[i - 1] > adjusted[i]) {
            adjusted[i - 1] = adjusted[i];
        }
    }

    for (size_t i = 0; i < n; ++i) {
        std::get<1>(data[i]) = std::min(1.0, adjusted[i]);
    }

    std::sort(data.begin(), data.end(), [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
    });
}

// Main Function
void add_BH_adjusted_column(
    const std::string& input_file, 
    const std::string& output_file_significant,
    const std::string& phenotype_type) {

    std::ifstream infile(input_file);
    std::string col;

    // First pass: Collect p-values
    std::vector<std::tuple<double, double, size_t>> pvalues;
    std::string line;
    size_t line_index = 0;
    size_t adjusted_col_index = phenotype_type == "binary" ? 6 : 5; 

    // Read the header line
    std::string header_line;
    std::getline(infile, header_line);
    std::stringstream header_ss(header_line);
    std::vector<std::string> headers;
    while (std::getline(header_ss, col, '\t')) {
        headers.push_back(col);
    }

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        double pval = 1.0;
        if (phenotype_type == "binary") {
            pval = mean_pvalue_from_strings(columns[4], columns[5]);
        } else if (phenotype_type == "quantitative") {
            pval = string_to_pvalue(columns[4]);
        }

        pvalues.emplace_back(pval, 1.0, line_index++);
    }
    infile.close();

    // Apply BH correction
    adjust_pvalues_BH(pvalues);

    // Second pass: rewrite with BH-adjusted values
    infile.open(input_file);
    std::ofstream outfile("temp_output.tsv");
    std::ofstream outfile_significant(output_file_significant);

    // Write headers
    outfile << header_line << '\n';
    outfile_significant << header_line << '\n';

    std::getline(infile, line); // Skip header again
    line_index = 0;

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> columns;

        while (std::getline(ss, token, '\t')) {
            columns.push_back(token);
        }

        double adjusted_p = std::get<1>(pvalues[line_index]);
        std::string adj_str = set_precision(adjusted_p);
        columns[adjusted_col_index] = adj_str;

        // Write updated line
        for (size_t i = 0; i < columns.size(); ++i) {
            outfile << columns[i];
            if (i != columns.size() - 1) outfile << '\t';
        }
        outfile << '\n';

        if (adjusted_p < 1e-5) {
            for (size_t i = 0; i < columns.size(); ++i) {
                outfile_significant << columns[i];
                if (i != columns.size() - 1) outfile_significant << '\t';
            }
            outfile_significant << '\n';
        }

        ++line_index;
    }

    infile.close();
    outfile.close();
    outfile_significant.close();

    // Replace original file
    std::remove(input_file.c_str());
    std::rename("temp_output.tsv", input_file.c_str());

}