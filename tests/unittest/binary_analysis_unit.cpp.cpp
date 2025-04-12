#include <catch2/catch_test_macros.hpp>
#include "../../src/binary_analysis.hpp"

TEST_CASE("Chi-square & Fisher test function", "[chi2Test]") {
    SECTION("Valid chi-square test & not valid Fisher test calculation") {
        std::vector<std::vector<size_t>> observed = {
            {10, 20, 30},
            {20, 30, 40}
        };
        REQUIRE(chi2Test(observed) != "");
        REQUIRE(fastFishersExactTest(observed) == "NA");
    }

    SECTION("Chi-square fail & Fisher test valid (zero row)") {
        std::vector<std::vector<size_t>> zero_row_matrix = {
            {0, 0},
            {4, 5}
        };
        REQUIRE(chi2Test(zero_row_matrix) == "1.0000");
        REQUIRE(fastFishersExactTest(zero_row_matrix) != "");
    }
}
