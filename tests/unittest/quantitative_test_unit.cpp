#include <catch.hpp>

#include "../../src/quantitative_table.hpp"
#include "../../src/stats_test.hpp"
#include "../../src/arg_parser.hpp"  // Qtl_data

// TEST_CASE("Linear Regression Test without cov", "[linear_regression]") {
//     SECTION("Linear Regression 1 - Perfect Linear Relationship") {

//         std::vector<std::vector<double>> df = {
//             {0, 1},
//             {1, 0},
//             {0, 0.5}
//         };

//         std::vector<double> quantitative_phenotype = {2.0, 4.0, 6.0};
//         std::vector<std::vector<double>> covar;  // No covariates

//         std::string se, beta, p_value, r2;

//         linear_regression(df, quantitative_phenotype, covar, p_value, beta, se, r2);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(se == "NA");
//         REQUIRE(beta == "-8.0000");
//         REQUIRE(p_value == "NA");
//         REQUIRE(r2 == "1.0000");
//     }

//     SECTION("Linear Regression 2 - Moderate") {

//         std::vector<std::vector<double>> df = {
//             {0.5, 0, 0.5},
//             {0, 0.5, 0.5},
//             {1, 0, 0},
//             {0, 1, 0},
//             {0, 0.5, 0}
//         };

//         std::vector<double> quantitative_phenotype = {10.5, 13.0, 15.8, 19.7, 21.5};
//         std::vector<std::vector<double>> covar;
//         std::string se, beta, p_value, r2;

//         linear_regression(df, quantitative_phenotype, covar, p_value, beta, se, r2);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(std::stod(se) == 0.880);
//         REQUIRE(std::stod(beta) == -17.4400);
//         REQUIRE(std::stod(p_value) == 0.0320);
//         REQUIRE(std::stod(r2) == 0.999);
//     }

//     SECTION("Linear Regression 3 - Weaker Correlation") {

//         std::vector<std::vector<double>> df = {
//             {1, 0, 0},
//             {1, 0, 0},
//             {1, 0, 0},
//             {1, 0, 0},
//             {1, 0, 0},
//             {1, 0, 0},
//             {1, 0, 0},
//             {0, 1, 0},
//             {0, 0, 0.5},

//         };

//         std::vector<double> quantitative_phenotype = {4.5, 7.0, 9.2, 10.9, 13.0, 14.0, 11.0, 15.0, 16.0};
//         std::vector<std::vector<double>> covar;  // No covariates
//         std::string se, beta, p_value, r2;
//         linear_regression(df, quantitative_phenotype, covar, p_value, beta, se, r2);

//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         REQUIRE(std::stod(se) == 3.033);
//         REQUIRE(std::stod(beta) == 6.5878);
//         REQUIRE(std::stod(p_value) == 0.0730);
//         REQUIRE(std::stod(r2) == 0.4210);
//     }
// }
