#include <catch2/catch_test_macros.hpp>
#include "../../src/matrix.hpp"


using namespace stoat_vcf;

TEST_CASE("stoat_vcf::EdgeBySampleMatrix Constructor and Basic Properties", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix initializes correctly") {
        std::vector<string> sample_names;
        stoat_vcf::EdgeBySampleMatrix mat(sample_names, 4, 5);
        REQUIRE(mat.get_matrix().size() > 0);  // Ensure matrix is allocated
        REQUIRE_FALSE(mat(0, 0));  // Initially, all elements should be false
    }
}

TEST_CASE("stoat_vcf::EdgeBySampleMatrix Expansion", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix expands properly") {
        std::vector<string> sample_names;
        stoat_vcf::EdgeBySampleMatrix mat(sample_names, 4, 5);
        size_t original_size = mat.get_matrix().size();
        mat.expandMatrix();
        REQUIRE(mat.get_matrix().size() > original_size);
    }
}
TEST_CASE("stoat_vcf::EdgeBySampleMatrix Set and Access Elements", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix correctly sets and retrieves values") {
        std::vector<string> sample_names;
        stoat_vcf::EdgeBySampleMatrix mat(sample_names, 4, 5);
        REQUIRE_FALSE(mat(1, 3));  // Initially, should be false
        mat.set(1, 3);
        REQUIRE(mat(1, 3));  // Should now be true
    }
}
TEST_CASE("stoat_vcf::EdgeBySampleMatrix Shrink", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix correctly shrinks") {
        std::vector<string> sample_names;
        stoat_vcf::EdgeBySampleMatrix mat(sample_names, 10, 5);
        size_t original_size = mat.get_matrix().size();
        mat.shrink(5);  // Reduce row count
        REQUIRE(mat.get_matrix().size() < original_size);
    }
}

TEST_CASE("stoat_vcf::EdgeBySampleMatrix Maximum Element", "[stoat_vcf::EdgeBySampleMatrix]") {
    SECTION("stoat_vcf::EdgeBySampleMatrix tracks maximum element correctly") {
        std::vector<string> sample_names;
        stoat_vcf::EdgeBySampleMatrix mat(sample_names, 4, 5);
        REQUIRE(mat.getMaxElement() == 4);  // Initially zero
        mat.expandMatrix();
        REQUIRE(mat.getMaxElement() == 8);  // Should be updated
    }
}
