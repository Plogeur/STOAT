#define CATCH_CONFIG_MAIN unit_test::binary
#include "catch.hpp"
#include "binary.hpp"

TEST_CASE("Matrix Initialization", "[Matrix]") {
    SECTION("Default initialization") {
        Matrix mat(3, 3);
        REQUIRE(mat.getRows() == 3);
        REQUIRE(mat.get_matrix().size() == 9);  // 3x3 flattened matrix
    }

    SECTION("Check default values") {
        Matrix mat(2, 2);
        REQUIRE_FALSE(mat(0, 0));  // Expect all to be false initially
        REQUIRE_FALSE(mat(1, 1));
    }
}

TEST_CASE("Matrix Set and Access", "[Matrix]") {
    Matrix mat(3, 3);

    SECTION("Set and retrieve values") {
        mat.set(0, 1);
        mat.set(2, 2);
        REQUIRE(mat(0, 1));
        REQUIRE(mat(2, 2));
        REQUIRE_FALSE(mat(1, 1));  // Ensure no unintended modifications
    }
}

TEST_CASE("Matrix Row Header Manipulation", "[Matrix]") {
    Matrix mat(2, 2);
    std::unordered_map<std::string, size_t> rowHeader = {{"A", 0}, {"B", 1}};
    
    mat.set_row_header(rowHeader);
    REQUIRE(mat.get_row_header().size() == 2);
    REQUIRE(mat.get_row_header().at("A") == 0);
    REQUIRE(mat.get_row_header().at("B") == 1);
}

TEST_CASE("Matrix Expansion", "[Matrix]") {
    Matrix mat(2, 2);
    size_t oldSize = mat.get_matrix().size();

    mat.expandMatrix();
    REQUIRE(mat.get_matrix().size() > oldSize);
}

TEST_CASE("Matrix Shrinking", "[Matrix]") {
    Matrix mat(5, 5);
    mat.shrink(3);
    REQUIRE(mat.getRows() == 3);
}

TEST_CASE("Binary Statistical Tests", "[Binary]") {
    SECTION("Chi-Square Test with valid data") {
        Matrix mat(2, 2);
        mat.set(0, 0);
        mat.set(1, 1);
        REQUIRE_NOTHROW(chi2Test(mat));
    }

    SECTION("Fisher's Exact Test on 2x2") {
        Matrix mat(2, 2);
        mat.set(0, 0);
        mat.set(1, 1);
        REQUIRE_NOTHROW(fastFishersExactTest(mat));
    }
}

TEST_CASE("Utility Functions", "[Binary]") {
    SECTION("String Join Function") {
        std::vector<std::string> words = {"one", "two", "three"};
        REQUIRE(join(words, "-") == "one-two-three");
    }
}
