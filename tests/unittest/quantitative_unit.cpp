#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"
#include "binary.hpp"

TEST_CASE("Mean Calculation", "[Statistics]") {
    SECTION("Computes correct mean for a sample vector") {
        std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0};
        REQUIRE(mean(values) == Approx(3.0));
    }
}

TEST_CASE("Variance Calculation", "[Statistics]") {
    SECTION("Computes correct variance for a sample vector") {
        std::vector<double> values = {2.0, 4.0, 6.0, 8.0, 10.0};
        REQUIRE(variance(values) == Approx(10.0));
    }
}

TEST_CASE("Covariance Calculation", "[Statistics]") {
    SECTION("Computes correct covariance for two sample vectors") {
        std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
        std::vector<double> y = {2.0, 4.0, 6.0, 8.0, 10.0};
        REQUIRE(covariance(x, y) == Approx(10.0));
    }

    SECTION("Throws exception for mismatched vector sizes") {
        std::vector<double> x = {1.0, 2.0, 3.0};
        std::vector<double> y = {4.0, 5.0};
        REQUIRE_THROWS_AS(covariance(x, y), std::invalid_argument);
    }
}

TEST_CASE("Linear Regression", "[Regression]") {
    SECTION("Computes correct regression parameters") {
        std::unordered_map<std::string, std::vector<int>> df = {
            {"A", {1, 2, 3}},
            {"B", {2, 4, 6}}
        };
        std::unordered_map<std::string, double> phenotype = {
            {"A", 10.0},
            {"B", 20.0}
        };

        auto [se, beta, p_val] = linear_regression(df, phenotype);
        
        REQUIRE(beta == Approx(5.0));  // Expected slope
        REQUIRE(se > 0);
    }

    SECTION("Returns NA if data is insufficient") {
        std::unordered_map<std::string, std::vector<int>> df = {{"A", {1}}};
        std::unordered_map<std::string, double> phenotype = {{"A", 10.0}};

        auto [se, beta, p_val] = linear_regression(df, phenotype);

        REQUIRE(p_val == "NA");
    }
}

TEST_CASE("Quantitative Table Creation", "[Table]") {
    SECTION("Creates correct genotype table") {
        Matrix mat(4, 4);
        mat.set_row_header({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 3}});

        std::vector<std::string> samples = {"A", "B", "C", "D"};
        std::vector<std::string> headers = {"A-B", "C-D"};

        auto df = create_quantitative_table(samples, headers, mat);
        
        REQUIRE(df.size() == 4);
        REQUIRE(df["A"].size() == 2);
    }
}

