#include <catch2/catch_test_macros.hpp>
#include "../../src/binary_analysis.hpp"

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

TEST_CASE("Chi-square & Fisher test function", "[chi2Test]") {
    SECTION("Valid chi-square test & not valid Fisher test calculation") {
        std::vector<std::vector<int>> observed = {
            {10, 20, 30},
            {20, 30, 40}
        };
        REQUIRE(chi2Test(observed) != "NA");
        REQUIRE(fastFishersExactTest(observed) == "NA");
    }

    SECTION("Chi-square fail & Fisher test valid (zero row)") {
        std::vector<std::vector<int>> zero_row_matrix = {
            {0, 0},
            {4, 5}
        };
        REQUIRE(chi2Test(zero_row_matrix) == "NA");
        REQUIRE(fastFishersExactTest(zero_row_matrix) != "NA");
    }
}
