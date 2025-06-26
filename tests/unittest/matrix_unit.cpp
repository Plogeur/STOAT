#include <catch2/catch_test_macros.hpp>
#include "../../src/matrix.hpp"

TEST_CASE("EdgeBySampleMatrix Constructor and Basic Properties", "[EdgeBySampleMatrix]") {
    SECTION("EdgeBySampleMatrix initializes correctly") {
        EdgeBySampleMatrix mat(4, 5);
        REQUIRE(mat.get_matrix().size() > 0);  // Ensure matrix is allocated
        REQUIRE_FALSE(mat(0, 0));  // Initially, all elements should be false
    }
}

TEST_CASE("EdgeBySampleMatrix Expansion", "[EdgeBySampleMatrix]") {
    SECTION("EdgeBySampleMatrix expands properly") {
        EdgeBySampleMatrix mat(4, 5);
        size_t original_size = mat.get_matrix().size();
        mat.expandMatrix();
        REQUIRE(mat.get_matrix().size() > original_size);
    }
}
TEST_CASE("EdgeBySampleMatrix Set and Access Elements", "[EdgeBySampleMatrix]") {
    SECTION("EdgeBySampleMatrix correctly sets and retrieves values") {
        EdgeBySampleMatrix mat(4, 5);
        REQUIRE_FALSE(mat(1, 3));  // Initially, should be false
        mat.set(1, 3);
        REQUIRE(mat(1, 3));  // Should now be true
    }
}
TEST_CASE("EdgeBySampleMatrix Shrink", "[EdgeBySampleMatrix]") {
    SECTION("EdgeBySampleMatrix correctly shrinks") {
        EdgeBySampleMatrix mat(10, 5);
        size_t original_size = mat.get_matrix().size();
        mat.shrink(5);  // Reduce row count
        REQUIRE(mat.get_matrix().size() < original_size);
    }
}

TEST_CASE("EdgeBySampleMatrix Maximum Element", "[EdgeBySampleMatrix]") {
    SECTION("EdgeBySampleMatrix tracks maximum element correctly") {
        EdgeBySampleMatrix mat(4, 5);
        REQUIRE(mat.getMaxElement() == 4);  // Initially zero
        mat.expandMatrix();
        REQUIRE(mat.getMaxElement() == 8);  // Should be updated
    }
}