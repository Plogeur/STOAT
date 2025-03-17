#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "../../src/quantitative_analysis.hpp"
#include "../../src/matrix.hpp"

using Catch::Approx;

// TEST_CASE("Linear Regression Test", "[linear_regression]") {
//     SECTION("Régression linéaire simple") {

//         std::unordered_map<std::string, std::vector<int>> df = {
//             {"Sample1", {1, 10, 20}},
//             {"Sample2", {2, 15, 25}},
//             {"Sample3", {3, 30, 35}}
//         };

//         std::unordered_map<std::string, double> quantitative_phenotype = {
//             {"Sample1", 2.0},
//             {"Sample2", 4.0},
//             {"Sample3", 6.0}
//         };

//         auto [se, beta, p_value, r2] = linear_regression(df, quantitative_phenotype);

//         // Afficher les valeurs pour le débogage
//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         // Vérifier que les valeurs sont correctes
//         REQUIRE(se != "NA");
//         REQUIRE(beta == "2.000");  // La pente devrait être exactement 2 (Y = 2X)
//         REQUIRE(p_value != "NA");
//         REQUIRE(r2 == "1.000");  // R² devrait être exactement 1 pour une relation linéaire parfaite
//     }

//     SECTION("Régression linéaire imparfaite") {
//         // Note: La régression utilise uniquement la première valeur de chaque vecteur
//         // X = [1, 2, 3, 4] et Y = [2, 3.9, 6.1, 7.8]
//         // Cela donne une relation approximativement linéaire avec :
//         // - pente proche de 2
//         // - R² < 1 car les points ne sont pas parfaitement alignés
//         std::unordered_map<std::string, std::vector<int>> df = {
//             {"Sample1", {1, 10, 20}},
//             {"Sample2", {2, 15, 25}},
//             {"Sample3", {3, 30, 35}},
//             {"Sample4", {4, 40, 45}}
//         };

//         std::unordered_map<std::string, double> quantitative_phenotype = {
//             {"Sample1", 2.0},
//             {"Sample2", 3.9},
//             {"Sample3", 6.1},
//             {"Sample4", 7.8}
//         };

//         auto [se, beta, p_value, r2] = linear_regression(df, quantitative_phenotype);

//         // Afficher les valeurs pour le débogage
//         INFO("se = " << se);
//         INFO("beta = " << beta);
//         INFO("p_value = " << p_value);
//         INFO("r2 = " << r2);

//         // Vérifier que les valeurs sont correctes
//         REQUIRE(se != "NA");  // Il devrait y avoir une erreur standard non nulle
//         REQUIRE(std::stod(beta) == Approx(1.93).margin(0.01));  // La pente devrait être proche de 1.93
//         REQUIRE(std::stod(p_value) < 0.05);  // La relation devrait être significative
//         REQUIRE(std::stod(r2) > 0.95);  // R² devrait être élevé mais pas égal à 1
//     }
// }

// TEST_CASE("Création de table quantitative", "[create_quantitative_table]") {
//     SECTION("Table simple") {
//         std::vector<std::string> list_samples = {"Sample1", "Sample2", "Sample3"};
//         std::vector<std::string> column_headers = {"Path1", "Path2"};

//         Matrix matrix(3, 2);  // 3 échantillons, 2 chemins
//         matrix.set(0, 0);  // Sample1, Path1
//         matrix.set(1, 1);  // Sample2, Path2

//         auto [table, size] = create_quantitative_table(list_samples, column_headers, matrix);

//         REQUIRE(table.size() == list_samples.size());
//         for (const auto& [sample, values] : table) {
//             REQUIRE(values.size() == column_headers.size());
//         }
//     }

//     SECTION("Table vide") {
//         std::vector<std::string> list_samples;
//         std::vector<std::string> column_headers;
//         Matrix matrix(0, 0);

//         auto [table, size] = create_quantitative_table(list_samples, column_headers, matrix);
//         REQUIRE(table.empty());
//     }
// }

// TEST_CASE("Test de set_precision", "[set_precision]") {
//     SECTION("Valeurs normales") {
//         REQUIRE(set_precision(3.14159) == "3.142");
//         REQUIRE(set_precision(0.0) == "0.000");
//         REQUIRE(set_precision(100.0) == "100.000");
//     }

//     SECTION("Valeurs extrêmes") {
//         REQUIRE(set_precision(1e-10) == "0.000");
//         REQUIRE(set_precision(1e10) == "10000000000.000");
//     }
// }