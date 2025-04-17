#include <catch2/catch_test_macros.hpp>
#include "../../src/binary_analysis.hpp"

TEST_CASE("Chi-square & Fisher test function", "[chi2Test]") {
    SECTION("Valid chi-square test & not valid Fisher test calculation") {
        std::vector<size_t> g0 = {10, 20, 30};
        std::vector<size_t> g1 = {20, 30, 40};
        REQUIRE(chi2Test(g0, g1) != "");
        REQUIRE(fastFishersExactTest(g0, g1) == "NA");
    }

    SECTION("Chi-square fail & Fisher test valid (zero row)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {20, 30};
        REQUIRE(chi2Test(g0, g1) == "1.0000");
        REQUIRE(fastFishersExactTest(g0, g1) != "");
    }
}
