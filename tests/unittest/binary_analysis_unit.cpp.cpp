#include <catch2/catch_test_macros.hpp>
#include "../../src/binary_analysis.hpp"

TEST_CASE("Chi-square & Fisher test function", "[chi2_2xN]") {
    SECTION("Valid chi-square test & valid Fisher test calculation") {
        std::vector<size_t> g0 = {10, 20};
        std::vector<size_t> g1 = {20, 10};
        REQUIRE(chi2_2x2(g0, g1) == "0.0201");
        REQUIRE(fastFishersExactTest(g0, g1) == "0.0194");
    }

    SECTION("Chi-square & Fisher test (significatif)") {
        std::vector<size_t> g0 = {30, 5};
        std::vector<size_t> g1 = {2, 25};
        REQUIRE(chi2_2x2(g0, g1) == "4.5938e-09");
        REQUIRE(fastFishersExactTest(g0, g1) == "3.5379e-10");
    }

    SECTION("Chi-square fail (N row)") {
        std::vector<size_t> g0 = {10, 15, 5};
        std::vector<size_t> g1 = {20, 10, 10};
        REQUIRE(chi2_2xN(g0, g1) == "0.0970");
    }

    SECTION("Chi-square fail (N row significatif)") {
        std::vector<size_t> g0 = {5, 10, 15, 20};
        std::vector<size_t> g1 = {20, 15, 10, 5};
        REQUIRE(chi2_2xN(g0, g1) == "0.0002");
    }

    SECTION("Chi-square fail (N row 1.0000)") {
        std::vector<size_t> g0 = {10, 10, 10, 10, 10};
        std::vector<size_t> g1 = {10, 10, 10, 10, 10};
        REQUIRE(chi2_2xN(g0, g1) == "1.0000");
    }

    SECTION("Chi-square fail & Fisher test fail (full zero row)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {0, 0};
        REQUIRE(chi2_2x2(g0, g1) == "NA");
        REQUIRE(fastFishersExactTest(g0, g1) == "NA");
    }

    SECTION("Chi-square fail (zero row)") {
        std::vector<size_t> g0 = {0, 0, 0};
        std::vector<size_t> g1 = {10, 20, 30};
        REQUIRE(chi2_2xN(g0, g1) == "NA");
    }

    SECTION("Chi-square fail & (zero column)") {
        std::vector<size_t> g0 = {0, 10, 5};
        std::vector<size_t> g1 = {0, 20, 15};
        REQUIRE(chi2_2xN(g0, g1) == "NA");
    }

    SECTION("Chi-square fail & Fisher test valid (zero row + column)") {
        std::vector<size_t> g0 = {0, 0};
        std::vector<size_t> g1 = {0, 1};
        REQUIRE(chi2_2x2(g0, g1) == "NA");
        REQUIRE(fastFishersExactTest(g0, g1) == "NA");
    }

    SECTION("Chi-square & Fisher test (zero / zero)") {
        std::vector<size_t> g0 = {1, 0};
        std::vector<size_t> g1 = {0, 1};
        REQUIRE(chi2_2x2(g0, g1) == "1.0000");
        REQUIRE(fastFishersExactTest(g0, g1) == "1.0000");
    }

    SECTION("Chi-square & Fisher test (strange but correct)") {
        std::vector<size_t> g0 = {79, 18};
        std::vector<size_t> g1 = {96, 23};
        REQUIRE(chi2_2x2(g0, g1) == "1.0000");
        REQUIRE(fastFishersExactTest(g0, g1) == "1.0000");
    }

    SECTION("Chi-square & Fisher test (very significative)") {
        std::vector<size_t> g0 = {122, 78};
        std::vector<size_t> g1 = {27, 173};
        REQUIRE(chi2_2x2(g0, g1) == "2.4445e-22");
        REQUIRE(fastFishersExactTest(g0, g1) == "1.4799e-23");
    }
}
