#include <catch.hpp>
#include "../../src/matrix.hpp"

using namespace stoat;

class TestEdgeBySampleMatrix : stoat::EdgeBySampleMatrix {
    public:
    TestEdgeBySampleMatrix(const std::vector<std::string>& sampleNames, size_t rows, size_t cols) : EdgeBySampleMatrix(sampleNames, rows, cols) {}

    using stoat::EdgeBySampleMatrix::matrix_1D;;
    using stoat::EdgeBySampleMatrix::operator();
    using stoat::EdgeBySampleMatrix::getMaxElement;
    using stoat::EdgeBySampleMatrix::expandMatrix;
    using stoat::EdgeBySampleMatrix::shrink;
    using stoat::EdgeBySampleMatrix::set;
};

TEST_CASE("stoat::EdgeBySampleMatrix Constructor and Basic Properties", "[stoat::EdgeBySampleMatrix]") {
    SECTION("stoat::EdgeBySampleMatrix initializes correctly") {
        std::vector<string> sample_names;
        TestEdgeBySampleMatrix mat(sample_names, 4, 5);
        REQUIRE(mat.matrix_1D.size() > 0);  // Ensure matrix is allocated
        REQUIRE_FALSE(mat(0, 0));  // Initially, all elements should be false
    }
}

TEST_CASE("stoat::EdgeBySampleMatrix Expansion", "[stoat::EdgeBySampleMatrix]") {
    SECTION("stoat::EdgeBySampleMatrix expands properly") {
        std::vector<string> sample_names;
        TestEdgeBySampleMatrix mat(sample_names, 4, 5);
        size_t original_size = mat.matrix_1D.size();
        mat.expandMatrix();
        REQUIRE(mat.matrix_1D.size() > original_size);
    }
}
TEST_CASE("stoat::EdgeBySampleMatrix Set and Access Elements", "[stoat::EdgeBySampleMatrix]") {
    SECTION("stoat::EdgeBySampleMatrix correctly sets and retrieves values") {
        std::vector<string> sample_names;
        TestEdgeBySampleMatrix mat(sample_names, 4, 5);
        REQUIRE_FALSE(mat(1, 3));  // Initially, should be false
        mat.set(1, 3);
        REQUIRE(mat(1, 3));  // Should now be true
    }
}
// TODO: Make this shrink to a specific size
TEST_CASE("stoat::EdgeBySampleMatrix Shrink", "[stoat::EdgeBySampleMatrix]") {
    SECTION("stoat::EdgeBySampleMatrix correctly shrinks") {
        std::vector<string> sample_names;
        TestEdgeBySampleMatrix mat(sample_names, 10, 5);
        size_t original_size = mat.matrix_1D.size();
        mat.shrink();  // Reduce row count
        REQUIRE(mat.matrix_1D.size() < original_size);
    }
}

TEST_CASE("stoat::EdgeBySampleMatrix Maximum Element", "[stoat::EdgeBySampleMatrix]") {
    SECTION("stoat::EdgeBySampleMatrix tracks maximum element correctly") {
        std::vector<string> sample_names;
        TestEdgeBySampleMatrix mat(sample_names, 4, 5);
        REQUIRE(mat.getMaxElement() == 4);  // Initially zero
        mat.expandMatrix();
        REQUIRE(mat.getMaxElement() == 8);  // Should be updated
    }
}
