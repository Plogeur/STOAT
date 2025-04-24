#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>

int main() {
    // Define the chi-squared distribution with 1 degree of freedom
    boost::math::chi_squared_distribution<long double> dist(1);

    // Loop over chi2_stat values from 20 to 50
    for (long double chi2_stat = 50.0; chi2_stat <= 300.0; chi2_stat += 1.0) {
        // Compute the p-value using the chi-squared CDF
        long double p_value = 1.0 - boost::math::cdf(dist, chi2_stat);

        // Print the chi2_stat and the corresponding p-value
        std::cout << "chi2_stat = " << chi2_stat << ", p_value = " << p_value << std::endl;
    }

    return 0;
}
