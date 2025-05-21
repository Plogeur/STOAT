#include <catch2/catch_test_macros.hpp>

#include "../../src/quantitative_analysis.hpp"
#include "../../src/matrix.hpp"

// TEST_CASE("create_quantitative_table basic behavior", "[quantitative]") {
//     SECTION("2 column simple case") {

//         size_t num_samples = 4;
//         std::vector<std::string> headers = {"path1", "path2"};
//         std::vector<double> phenotype = {-2.4234, -1.3242, 0.3214, 2.3248};

//         // Create a simple matrix (4 samples / 2 haplotypes)
//         Matrix matrix(4, 2);

//         // Simulate that paths match at certain indices
//         matrix.set(0, 1);
//         matrix.set(0, 1);
//         matrix.set(1, 0);
//         matrix.set(1, 0);
//         matrix.set(2, 1);
//         matrix.set(2, 1);
//         matrix.set(3, 1);
//         matrix.set(3, 0);

//         auto [genotypes, phenotype_filtered, allele_count] =
//             create_quantitative_table(num_samples, headers, phenotype, matrix);

//         REQUIRE(genotypes.size() == 4);
//         REQUIRE(genotypes[0].size() == 2);
//         REQUIRE(phenotype_filtered.size() == 3);
//         REQUIRE(allele_count == 5);
//         REQUIRE(phenotype_filtered[0] == -2.4234);
//         REQUIRE(phenotype_filtered[1] == 0.3214); // remove the second pheno without allele passing throught
//         REQUIRE(phenotype_filtered[2] == 2.3248);
//     }

//     SECTION("2 column case 1th col full 0") {
//         size_t num_samples = 4;
//         std::vector<std::string> headers = {"path1", "path2"};
//         std::vector<double> phenotype = {-2.4234, -1.3242, 0.3214, 2.3248};

//         // Create a simple matrix (4 samples / 2 haplotypes)
//         Matrix matrix(4, 2);

//         // Simulate that paths match at certain indices
//         matrix.set(0, 0);
//         matrix.set(0, 1);
//         matrix.set(1, 0);
//         matrix.set(1, 1);
//         matrix.set(2, 0);
//         matrix.set(2, 1);
//         matrix.set(3, 0);
//         matrix.set(3, 0);

//         auto [genotypes, phenotype_filtered, allele_count] =
//             create_quantitative_table(num_samples, headers, phenotype, matrix);

//         REQUIRE(genotypes.size() == 4);
//         REQUIRE(genotypes[0].size() == 1); // remove the 2 column composed of full 0
//         REQUIRE(phenotype_filtered.size() == 3);
//         REQUIRE(allele_count == 3);
//         REQUIRE(phenotype_filtered[0] == -2.4234);
//         REQUIRE(phenotype_filtered[1] == -1.3242);
//         REQUIRE(phenotype_filtered[2] == 0.3214); // remove the third pheno without allele passing throught
//     }

//     SECTION("3 column case") {
//         size_t num_samples = 4;
//         std::vector<std::string> headers = {"path1", "path2"};
//         std::vector<double> phenotype = {-2.4234, -1.3242, 0.3214, 2.3248};

//         // Create a simple matrix (4 samples / 3 haplotypes)
//         Matrix matrix(2, 8);

//         // Simulate that paths match at certain indices
//         matrix.set(0, 1);
//         matrix.set(0, 1);
//         matrix.set(0, 1);
//         matrix.set(0, 1);
//         matrix.set(1, 1);
//         matrix.set(1, 1);
//         matrix.set(1, 1);
//         matrix.set(1, 1);

//         auto [genotypes, phenotype_filtered, allele_count] =
//             create_quantitative_table(num_samples, headers, phenotype, matrix);

//         REQUIRE(genotypes.size() == 4);
//         REQUIRE(genotypes[0].size() == 2);
//         REQUIRE(phenotype_filtered.size() == 4);
//         REQUIRE(allele_count == 3);
//     }

//     SECTION("3 column case last col full 0") {

//     }
// }

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
// }

// TEST_CASE("Test de set_precision", "[set_precision]") {
//     SECTION("Valeurs normales") {
//         REQUIRE(set_precision(3.14159) == "3.1416");
//         REQUIRE(set_precision(0.0) == "0.0000e+00");
//         REQUIRE(set_precision(100.0) == "100.0000");
//     }

//     SECTION("Valeurs extrêmes") {
//         REQUIRE(set_precision(1e-10) == "1.0000e-10");
//         REQUIRE(set_precision(1e10) == "10000000000.0000");
//     }
// }