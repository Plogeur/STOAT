#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

// Holm-Bonferroni correction in C++
std::vector<double> stoat_vcf::adjusted_holm(const std::vector<double>& p_values) {
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

int main() {
    std::vector<double> p_values = {0.053257545459073324, 0.16349746926294687};
    std::vector<double> adjusted = adjusted_holm(p_values);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Original p-values:\n";
    for (double p : p_values) std::cout << p << " ";
    std::cout << "\n\nHolm-adjusted p-values:\n";
    for (double p : adjusted) std::cout << p << " ";
    std::cout << std::endl;

    return 0;
}
