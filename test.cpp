#include <iostream>
#include <chrono>
#include <cmath>

// ============================= First approach =============================
// Fisher's Exact Test for 2x2 contingency table
#ifndef DBL_MAX
#  define DBL_MAX 1.7976931348623157e308
#endif

#ifdef __cplusplus
#  define K_CAST(type, val) (const_cast<type>(val))
#  define R_CAST(type, val) (reinterpret_cast<type>(val))
#  define S_CAST(type, val) (static_cast<type>(val))
#endif

// 2^{-40} for now, since 2^{-44} was too small on real data
static const double kExactTestEpsilon2 = 0.0000000000009094947017729282379150390625;
static const double kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;

int main() {
    // Example 2x2 table values
    size_t m11 = 1, m12 = 9, m21 = 11, m22 = 5;

    // Run the first approach 1000 times and time it
    auto start1 = std::chrono::high_resolution_clock::now();
    double p_value1 = 0;
    for (int i = 0; i < 10000000; ++i) {
        p_value1 = FisherExact2x2P(m11, m12, m21, m22);
    }
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start1;
    std::cout << "First approach (1000000 runs) p-value: " << p_value1 << std::endl;
    std::cout << "First approach (1000000 runs) execution time: " << duration1.count() << " seconds" << std::endl;

    // Run the second approach 1000 times and time it
    auto start2 = std::chrono::high_resolution_clock::now();
    double p_value2 = 0;
    for (int i = 0; i < 10000000; ++i) {
        p_value2 = fastFishersExactTest(m11, m12, m21, m22);
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration2 = end2 - start2;
    std::cout << "Second approach (1000000 runs) p-value: " << p_value2 << std::endl;
    std::cout << "Second approach (1000000 runs) execution time: " << duration2.count() << " seconds" << std::endl;

    return 0;
}