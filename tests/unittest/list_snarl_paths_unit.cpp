#include <catch2/catch_test_macros.hpp>
#include "../../src/list_snarl_paths.hpp"

TEST_CASE("Test de la classe Path", "[Path]") {
    SECTION("Construction et ajout de nœuds") {
        Path path;
        path.addNode("1", '>');
        path.addNode("2", '<');
        path.addNode("3", '>');

        REQUIRE(path.size() == 3);
        REQUIRE(path.print() == ">1<2>3");
        REQUIRE(path.nreversed() == 1);
    }

    SECTION("Test de flip") {
        Path path;
        path.addNode("1", '>');
        path.addNode("2", '<');
        path.addNode("3", '>');

        path.flip();
        REQUIRE(path.print() == "<3>2<1");
        REQUIRE(path.nreversed() == 2);
    }

    SECTION("Chemin vide") {
        Path path;
        REQUIRE(path.size() == 0);
        REQUIRE(path.print() == "");
        REQUIRE(path.nreversed() == 0);
    }
}

TEST_CASE("Test de calcul_pos_type_variant", "[calcul_pos_type_variant]") {

    SECTION("Liste simple") {
        std::vector<std::tuple<std::string, size_t, size_t, size_t>> list_paths = {
            {"T", 0, 0, 3},  // SNP: "T", length 1
            {"TT", 0, 0, 3}, // INS: "TT", length 2
            {"", 0, 0, 2}    // DEL: "", length 0
        };
        auto [types, padding] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "T");  // SNP
        REQUIRE(types[1] == "2");  // INS
        REQUIRE(types[2] == "0");  // DEL
        REQUIRE(padding == 0);
    }

    SECTION("Cas SNP uniquement") {
        std::vector<std::tuple<std::string, size_t, size_t, size_t>> list_paths = {
            {"T", 0, 0, 3},
            {"G", 0, 0, 3}
        };
        auto [types, padding] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 2);
        REQUIRE(types[0] == "T");
        REQUIRE(types[1] == "G");
        REQUIRE(padding == 1);  // All SNPs → padding = 1
    }

    SECTION("Cas complexe") {
        std::vector<std::tuple<std::string, size_t, size_t, size_t>> list_paths = {
            {"_", 10, 20, 10},   // Complex path > 3
            {"TTTT", 0, 0, 3},   // Insertion
            {"", 0, 0, 2}        // Deletion
        };
        auto [types, padding] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "CPX:10/20");
        REQUIRE(types[1] == "4");
        REQUIRE(types[2] == "0");
        REQUIRE(padding == 0);
    }
}
