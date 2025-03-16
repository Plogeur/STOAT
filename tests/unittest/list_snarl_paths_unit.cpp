#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "list_snarl_paths.hpp"

using Catch::Approx;

TEST_CASE("Test de la classe Path", "[Path]") {
    SECTION("Construction et ajout de nœuds") {
        Path path;
        path.addNode("1", '+');
        path.addNode("2", '-');
        path.addNode("3", '+');
        
        REQUIRE(path.size() == 3);
        REQUIRE(path.print() == "1>2<3");
        REQUIRE(path.nreversed() == 1);
    }

    SECTION("Test de flip") {
        Path path;
        path.addNode("1", '+');
        path.addNode("2", '-');
        path.addNode("3", '+');
        
        path.flip();
        REQUIRE(path.print() == "3>2>1");
        REQUIRE(path.nreversed() == 0);
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
        auto [types, size] = calcul_pos_type_variant(empty_list);
        REQUIRE(types.empty());
        REQUIRE(size == 0);
    }

    SECTION("Liste simple") {
        std::vector<std::vector<std::string>> list_paths = {
            {"A", "T", "G"},  // SNP
            {"A", "TT", "G"}, // INS
            {"A", "G"}        // DEL
        };
        auto [types, size] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "T");
        REQUIRE(types[1] == "TT");
        REQUIRE(types[2] == "DEL");
        REQUIRE(size == 3);
    }

    SECTION("Cas SNP uniquement") {
        std::vector<std::vector<std::string>> list_paths = {
            {"A", "T", "G"},
            {"C", "G", "T"}
        };
        auto [types, size] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 2);
        REQUIRE(types[0] == "T");
        REQUIRE(types[1] == "G");
        REQUIRE(size == 2);
    }

    SECTION("Cas complexe") {
        std::vector<std::vector<std::string>> list_paths = {
            {"A", "_", "G"},  // Complex
            {"A", "TTTT", "G"}, // INS
            {"A", "G"}        // DEL
        };
        auto [types, size] = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "CPX");
        REQUIRE(types[1] == "INS");
        REQUIRE(types[2] == "DEL");
        REQUIRE(size == 3);
    }
}

TEST_CASE("Test de check_threshold", "[check_threshold]") {
    SECTION("Valeurs valides") {
        REQUIRE_NOTHROW(check_threshold(0.0));
        REQUIRE_NOTHROW(check_threshold(0.5));
        REQUIRE_NOTHROW(check_threshold(1.0));
    }

    SECTION("Valeurs invalides") {
        REQUIRE_THROWS_AS(check_threshold(-0.1), std::invalid_argument);
        REQUIRE_THROWS_AS(check_threshold(1.1), std::invalid_argument);
    }
}

TEST_CASE("Test de find_snarl_id", "[find_snarl_id]") {
    SECTION("Documentation du test") {
        WARN("Les tests de find_snarl_id nécessitent une configuration plus complexe avec des données de test réelles");
    }
}

TEST_CASE("Test de follow_edges", "[follow_edges]") {
    SECTION("Documentation du test") {
        WARN("Les tests de follow_edges nécessitent une configuration plus complexe avec des données de test réelles");
    }
}

TEST_CASE("Test de fill_pretty_paths", "[fill_pretty_paths]") {
    SECTION("Documentation du test") {
        WARN("Les tests de fill_pretty_paths nécessitent une configuration plus complexe avec des données de test réelles");
    }
}