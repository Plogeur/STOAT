#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "matrix.hpp"

using Catch::Approx;

TEST_CASE("Matrix Constructor and Basic Properties", "[Matrix]") {
    SECTION("Matrix initializes correctly") {
        Matrix mat(4, 5);
        REQUIRE(mat.get_matrix().size() > 0);  // Ensure matrix is allocated
        REQUIRE_FALSE(mat(0, 0));  // Initially, all elements should be false
    }
}

TEST_CASE("Matrix Getters and Setters", "[Matrix]") {
    SECTION("Row header is set and retrieved correctly") {
        Matrix mat(4, 5);
        std::unordered_map<std::string, size_t> headers = {{"A", 0}, {"B", 1}, {"C", 2}, {"D", 3}};
        mat.set_row_header(headers);
        REQUIRE(mat.get_row_header().size() == 4);
        REQUIRE(mat.get_row_header().at("A") == 0);
    }
}

TEST_CASE("Matrix Expansion", "[Matrix]") {
    SECTION("Matrix expands properly") {
        Matrix mat(4, 5);
        size_t original_size = mat.get_matrix().size();
        mat.expandMatrix();
        REQUIRE(mat.get_matrix().size() > original_size);
    }
}

TEST_CASE("Matrix Set and Access Elements", "[Matrix]") {
    SECTION("Matrix correctly sets and retrieves values") {
        Matrix mat(4, 5);
        REQUIRE_FALSE(mat(1, 3));  // Initially, should be false
        mat.set(1, 3);
        REQUIRE(mat(1, 3));  // Should now be true
    }
}

TEST_CASE("Matrix Shrink", "[Matrix]") {
    SECTION("Matrix correctly shrinks") {
        Matrix mat(10, 5);
        size_t original_size = mat.get_matrix().size();
        mat.shrink(5);  // Reduce row count
        REQUIRE(mat.get_matrix().size() < original_size);
    }
}

TEST_CASE("Matrix Maximum Element", "[Matrix]") {
    SECTION("Matrix tracks maximum element correctly") {
        Matrix mat(4, 5);
        REQUIRE(mat.getMaxElement() == 0);  // Initially zero
        mat.set(1, 3);  // Set an element
        REQUIRE(mat.getMaxElement() > 0);  // Should be updated
    }
}
