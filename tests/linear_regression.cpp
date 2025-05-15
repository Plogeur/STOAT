#include <iostream>
#include <vector>
#include <algorithm>
#include "dataanalysis.h"  // ALGLIB header

// --- Holm-Bonferroni Correction ---
std::vector<double> adjusted_holm(const std::vector<double>& p_values) {
    int m = p_values.size();
    std::vector<std::pair<double, int>> indexed;
    for (int i = 0; i < m; ++i) {
        indexed.emplace_back(p_values[i], i);
    }

    std::sort(indexed.begin(), indexed.end());

    std::vector<double> adjusted(m);
    double prev = 0.0;
    for (int i = 0; i < m; ++i) {
        double raw = (m - i) * indexed[i].first;
        raw = std::min(raw, 1.0);
        adjusted[i] = std::max(prev, raw);
        prev = adjusted[i];
    }

    std::vector<double> reordered(m);
    for (int i = 0; i < m; ++i) {
        reordered[indexed[i].second] = adjusted[i];
    }

    return reordered;
}

// Convert std::vector<std::vector<double>> to alglib::real_2d_array
alglib::real_2d_array vector_to_alglib_matrix(const std::vector<std::vector<double>>& data) {
    int rows = data.size();
    int cols = data[0].size();
    alglib::real_2d_array arr;
    arr.setlength(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            arr[i][j] = data[i][j];
    return arr;
}

// Convert std::vector<double> to alglib::real_1d_array
alglib::real_1d_array vector_to_alglib_array(const std::vector<double>& data) {
    alglib::real_1d_array arr;
    arr.setlength(data.size());
    for (int i = 0; i < (int)data.size(); ++i)
        arr[i] = data[i];
    return arr;
}

void linear_regression(const std::vector<std::vector<double>>& features,
                       const std::vector<double>& phenotype) {

    int n_samples = features.size();
    int n_features = features[0].size();

    alglib::real_2d_array X = vector_to_alglib_matrix(features);
    alglib::real_1d_array y = vector_to_alglib_array(phenotype);

    alglib::real_1d_array c;   // regression coefficients
    alglib::lsfitreport rep;   // fit report (statistics)

    // Perform linear regression with intercept
    alglib::lsfitlinear(X, y, n_samples, n_features, c, rep, true);

    std::cout << "Linear Regression Results:\n";
    std::cout << "R²: " << rep.rsq << "\n\n";

    std::cout << "Feature\tBeta\t\tStdErr\t\tT-stat\t\tP-value\t\tP-adj\n";

    // Collect p-values for Holm correction (skip intercept)
    std::vector<double> p_values;
    for (int i = 1; i < c.length(); ++i) {
        p_values.push_back(rep.pvalues[i]);
    }
    std::vector<double> p_adj = adjusted_holm(p_values);

    // Print intercept (no adjustment)
    std::cout << "Intercept\t" << c[0] << "\t" << rep.stderr[0] << "\t" << rep.tvalues[0] << "\t" << rep.pvalues[0] << "\t" << "-" << "\n";

    // Print features with adjusted p-values
    for (int i = 1; i < c.length(); ++i) {
        std::cout << "X" << i << "\t\t" << c[i] << "\t" << rep.stderr[i] << "\t" << rep.tvalues[i] << "\t" << rep.pvalues[i] << "\t" << p_adj[i - 1] << "\n";
    }
}

int main() {
    // --- Synthetic Data ---
    std::vector<std::vector<double>> features = {
        {1,  5,  8},
        {2,  6,  7},
        {3,  5,  6},
        {4,  5,  5},
        {5,  6,  4},
        {6,  5,  3},
        {7,  4,  2},
        {8,  6,  1},
        {9,  5,  0},
        {10, 4, -1},
        {11, 5, -2},
        {12, 4, -3},
        {13, 5, -4},
        {14, 6, -5},
        {15, 5, -6},
        {16, 4, -7},
        {17, 6, -8},
        {18, 5, -9},
        {19, 4, -10},
        {20, 5, -11}
    };

    std::vector<double> phenotype = {
        2.1, 3.2, 4.1, 5.3, 6.0,
        7.2, 8.1, 9.4, 10.5, 11.7,
        12.9, 14.0, 14.8, 16.1, 17.4,
        18.2, 19.6, 20.5, 21.9, 22.7
    };

    linear_regression(features, phenotype);

    return 0;
}
