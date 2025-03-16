#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

using namespace std;

std::string set_precision(double value) {
    std::ostringstream oss;
    if (std::abs(value) < 0.0001) {
        return "0.000";
    } else {
        oss << std::fixed << std::setprecision(3) << value;
    }
    return oss.str();
}

// Linear regression function OLS
std::tuple<string, string, string, string> linear_regression(
    const std::unordered_map<std::string, std::vector<int>>& df,
    const std::unordered_map<std::string, double>& quantitative_phenotype) {
    
    try {
        // Check for missing data
        for (const auto& [sample, _] : df) {
            if (quantitative_phenotype.find(sample) == quantitative_phenotype.end()) {
                return {"NA", "NA", "NA", "NA"};
            }
        }
        for (const auto& [sample, _] : quantitative_phenotype) {
            if (df.find(sample) == df.end()) {
                return {"NA", "NA", "NA", "NA"};
            }
        }

        size_t num_samples = df.size();
        if (num_samples < 3) {  // Au moins 3 observations pour une régression
            return {"NA", "NA", "NA", "NA"};
        }

        // Vérifier que tous les vecteurs ont la même taille
        size_t vec_size = df.begin()->second.size();
        for (const auto& [_, vec] : df) {
            if (vec.size() != vec_size) {
                return {"NA", "NA", "NA", "NA"};
            }
        }

        // Utiliser Eigen pour la régression
        Eigen::MatrixXd X(num_samples, vec_size);
        Eigen::VectorXd y(num_samples);
        
        int row = 0;
        for (const auto& [sample, paths] : df) {
            y(row) = quantitative_phenotype.at(sample);
            for (size_t col = 0; col < paths.size(); ++col) {
                X(row, col) = static_cast<double>(paths[col]);
            }
            ++row;
        }

        // Calculer les coefficients de régression
        Eigen::VectorXd beta = (X.transpose() * X).ldlt().solve(X.transpose() * y);
        Eigen::VectorXd y_pred = X * beta;
        Eigen::VectorXd residuals = y - y_pred;

        // Calculer R²
        double rss = residuals.squaredNorm();
        double tss = (y.array() - y.mean()).matrix().squaredNorm();
        double r2 = 1.0 - (rss / tss);

        // Calculer l'erreur standard
        double mse = rss / (num_samples - vec_size);
        Eigen::MatrixXd cov_matrix = mse * (X.transpose() * X).inverse();
        double se = std::sqrt(cov_matrix(0,0));

        // Calculer la statistique F et la p-value
        double f_stat = (r2 / (vec_size - 1)) / ((1.0 - r2) / (num_samples - vec_size));
        double p_value = 1.0;
        if (f_stat > 0) {
            boost::math::fisher_f dist(vec_size - 1, num_samples - vec_size);
            p_value = boost::math::cdf(boost::math::complement(dist, f_stat));
        }

        return {
            set_precision(se),
            set_precision(beta(0)),
            set_precision(p_value),
            set_precision(r2)
        };
    } catch (...) {
        return {"NA", "NA", "NA", "NA"};
    }
}

// Function to create the quantitative table
std::pair<std::unordered_map<std::string, std::vector<int>>, size_t> create_quantitative_table(
    const std::vector<std::string>& list_samples, 
    const std::vector<std::string>& column_headers,
    Matrix& matrix) {

    // Check for empty inputs
    if (list_samples.empty() || column_headers.empty()) {
        return {std::unordered_map<std::string, std::vector<int>>(), 0};
    }

    // Retrieve row headers dictionary
    std::unordered_map<std::string, size_t> row_headers_dict = matrix.get_row_header();
    size_t allele_number = 0;
    size_t length_sample = list_samples.size();
    size_t length_column = column_headers.size();

    // Initialize a zero matrix for genotypes
    std::vector<std::vector<int>> genotypes(length_sample, std::vector<int>(length_column, 0));

    try {
        // Genotype paths
        for (size_t col_idx = 0; col_idx < length_column; ++col_idx) {
            const std::string& path_snarl = column_headers[col_idx];
            std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

            // Identify correct paths
            std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, 
                                                                matrix, length_sample * 2);

            for (auto idx : idx_srr_save) {
                if (idx >= 0 && idx/2 < static_cast<int>(length_sample)) {
                    size_t srr_idx = idx / 2;
                    genotypes[srr_idx][col_idx] += 1;
                    allele_number++;
                }
            }
        }

        // Create the final map
        std::unordered_map<std::string, std::vector<int>> df;
        for (size_t i = 0; i < list_samples.size(); ++i) {
            df[list_samples[i]] = genotypes[i];
        }
        
        return {df, allele_number};
    } catch (...) {
        return {std::unordered_map<std::string, std::vector<int>>(), 0};
    }
}
