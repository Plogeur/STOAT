#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "../../src/binary_analysis.hpp"

using Catch::Approx;

TEST_CASE("Check observed matrix validity", "[check_observed]") {
    SECTION("Valid matrix (no zero rows or columns)") {
        std::vector<std::vector<int>> valid_matrix = {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        };
        REQUIRE(check_observed(valid_matrix, 3, 3) == true);
    }

    SECTION("Matrix with a zero row") {
        std::vector<std::vector<int>> zero_row_matrix = {
            {0, 0, 0},
            {4, 5, 6},
            {7, 8, 9}
        };
        REQUIRE(check_observed(zero_row_matrix, 3, 3) == false);
    }

    SECTION("Matrix with a zero column") {
        std::vector<std::vector<int>> zero_column_matrix = {
            {1, 2, 0},
            {4, 5, 0},
            {7, 8, 0}
        };
        REQUIRE(check_observed(zero_column_matrix, 3, 3) == false);
    }

    SECTION("Empty matrix") {
        std::vector<std::vector<int>> empty_matrix = {};
        REQUIRE(check_observed(empty_matrix, 0, 0) == false);
    }
}

TEST_CASE("Chi-square test function", "[chi2Test]") {
    SECTION("Valid chi-square test calculation") {
        std::vector<std::vector<int>> observed = {
            {10, 20, 30},
            {20, 30, 40},
            {30, 40, 50}
        };
        REQUIRE(chi2Test(observed) != "NA");
    }

    SECTION("Chi-square test on an invalid matrix (zero row)") {
        std::vector<std::vector<int>> zero_row_matrix = {
            {0, 0, 0},
            {4, 5, 6},
            {7, 8, 9}
        };
        REQUIRE(chi2Test(zero_row_matrix) == "NA");
    }

    SECTION("Chi-square test on an invalid matrix (zero column)") {
        std::vector<std::vector<int>> zero_column_matrix = {
            {1, 2, 0},
            {4, 5, 0},
            {7, 8, 0}
        };
        REQUIRE(chi2Test(zero_column_matrix) == "NA");
    }
}

TEST_CASE("Log Factorial Initialization", "[initLogFacs]") {
    int n = 5;
    long double logFacs[n+1]; 

    initLogFacs(logFacs, n);

    SECTION("Base case: logFacs[0] should be 0") {
        REQUIRE(logFacs[0] == Approx(0.0));
    }

    SECTION("Computed log factorials") {
        REQUIRE(logFacs[1] == Approx(std::log(1.0)));
        REQUIRE(logFacs[2] == Approx(std::log(1.0) + std::log(2.0)));
        REQUIRE(logFacs[3] == Approx(std::log(1.0) + std::log(2.0) + std::log(3.0)));
        REQUIRE(logFacs[4] == Approx(std::log(1.0) + std::log(2.0) + std::log(3.0) + std::log(4.0)));
    }
}

TEST_CASE("Log Hypergeometric Probability Calculation", "[logHypergeometricProb]") {
    int n = 10;
    long double logFacs[n+1]; 
    initLogFacs(logFacs, n);

    SECTION("Basic case: small contingency table") {
        int a = 2, b = 3, c = 4, d = 1;
        long double result = logHypergeometricProb(logFacs, a, b, c, d);

        REQUIRE(std::isfinite(result)); // Ensure it's a valid number
    }

    SECTION("Edge case: all zeros (should return zero probability)") {
        int a = 0, b = 0, c = 0, d = 0;
        long double result = logHypergeometricProb(logFacs, a, b, c, d);

        REQUIRE(result == Approx(0.0));
    }
}

TEST_CASE("Fast Fisher's Exact Test", "[fastFishersExactTest]") {
    SECTION("Valid 2x2 contingency table") {
        std::vector<std::vector<int>> table = {
            {1, 9},
            {11, 3}
        };

        long double result = fastFishersExactTest(table);

        REQUIRE(result >= 0.0);
        REQUIRE(result <= 1.0);
    }

    SECTION("Invalid table (not 2x2)") {
        std::vector<std::vector<int>> invalid_table = {
            {1, 2, 3},
            {4, 5, 6}
        };

        REQUIRE(fastFishersExactTest(invalid_table) == Approx(-1.0));
    }

    SECTION("Edge case: all zeros (should return -1)") {
        std::vector<std::vector<int>> zero_table = {
            {0, 0},
            {0, 0}
        };

        REQUIRE(fastFishersExactTest(zero_table) == Approx(-1.0));
    }
}

TEST_CASE("Format group paths correctly", "[format_group_paths]") {
    SECTION("Empty input returns empty string") {
        std::vector<std::vector<int>> empty_df;
        REQUIRE(format_group_paths(empty_df) == "");
    }

    SECTION("Single row and column") {
        std::vector<std::vector<int>> df = {{1}};
        REQUIRE(format_group_paths(df) == "1");
    }

    SECTION("Multiple rows and columns") {
        std::vector<std::vector<int>> df = {
            {1, 2, 3},
            {4, 5, 6}
        };
        REQUIRE(format_group_paths(df) == "1:4,2:5,3:6");
    }
}

TEST_CASE("Binary Stat Test computes correctly", "[binary_stat_test]") {
    SECTION("Valid contingency table") {
        std::vector<std::vector<int>> table = {
            {1, 2},
            {3, 4}
        };
        std::vector<std::string> results = binary_stat_test(table);

        REQUIRE(results.size() == 8);
        REQUIRE(results[0] != "NA");  // Fisher's p-value should be a number
        REQUIRE(std::stod(results[2]) > 0); // Allele number
        REQUIRE(std::stod(results[3]) > 0); // Min row index
        REQUIRE(std::stod(results[4]) == 2); // Number of columns
        REQUIRE(std::stod(results[5]) >= 0); // Inter group
        REQUIRE(std::stod(results[6]) > 0); // Average should be valid
    }

    SECTION("Empty input should return default values") {
        std::vector<std::vector<int>> empty_df;
        std::vector<std::string> results = binary_stat_test(empty_df);

        REQUIRE(results.size() == 8);
        REQUIRE(results[0] == "NA");
        REQUIRE(results[2] == "0");
        REQUIRE(results[3] == std::to_string(std::numeric_limits<int>::max()));
        REQUIRE(results[4] == "0");
        REQUIRE(results[5] == "0");
        REQUIRE(results[6] == "0.000000");
        REQUIRE(results[7] == "");
    }
}
