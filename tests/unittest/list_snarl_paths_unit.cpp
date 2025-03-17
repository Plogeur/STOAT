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
    SECTION("Liste vide") {
        std::vector<std::vector<std::string>> empty_list;
        auto [types, padding] = calcul_pos_type_variant(empty_list);
        REQUIRE(types.empty());
    }

    SECTION("Liste simple") {
        std::vector<std::vector<std::string>> list_paths = {
            {"A", "T", "G"},  // SNP
            {"A", "TT", "G"}, // INS
            {"A", "G"}        // DEL
        };
        auto [types, padding] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "T");
        REQUIRE(types[1] == "TT");
        REQUIRE(types[2] == "DEL");
        REQUIRE(padding == 0);
    }

    SECTION("Cas SNP uniquement") {
        std::vector<std::vector<std::string>> list_paths = {
            {"A", "T", "G"},
            {"C", "G", "T"}
        };
        auto [types, padding] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 2);
        REQUIRE(types[0] == "T");
        REQUIRE(types[1] == "G");
        REQUIRE(padding == 1);
    }

    SECTION("Cas complexe") {
        std::vector<std::vector<std::string>> list_paths = {
            {"A", "_", "G"},  // Complex
            {"A", "TTTT", "G"}, // INS
            {"A", "G"}        // DEL
        };
        auto [types, padding] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "CPX");
        REQUIRE(types[1] == "INS");
        REQUIRE(types[2] == "DEL");
        REQUIRE(padding == 0);
    }
}
