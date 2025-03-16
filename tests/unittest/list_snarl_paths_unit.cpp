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
        vector<vector<string>> empty_list;
        auto [positions, size] = calcul_pos_type_variant(empty_list);
        REQUIRE(positions.empty());
        REQUIRE(size == 0);
    }

    SECTION("Liste simple") {
        vector<vector<string>> paths = {
            {"1", "2", "3"},
            {"1", "4", "3"},
            {"1", "5", "3"}
        };
        auto [positions, size] = calcul_pos_type_variant(paths);
        REQUIRE(positions.size() == 3);
        REQUIRE(size == 3);
    }
}

TEST_CASE("Test de check_threshold", "[check_threshold]") {
    SECTION("Proportion valide") {
        REQUIRE_NOTHROW(check_threshold(0.5));
        REQUIRE_NOTHROW(check_threshold(0.0));
        REQUIRE_NOTHROW(check_threshold(1.0));
    }

    SECTION("Proportion invalide") {
        REQUIRE_THROWS(check_threshold(-0.1));
        REQUIRE_THROWS(check_threshold(1.1));
    }
}

TEST_CASE("Test de find_snarl_id", "[find_snarl_id]") {
    // Note: Ces tests nécessitent une instance de SnarlDistanceIndex et net_handle_t
    // qui sont difficiles à créer dans un contexte de test unitaire.
    // Il faudrait créer des mocks ou utiliser des données de test réelles.
    SECTION("Documentation du test") {
        WARN("Les tests de find_snarl_id nécessitent une configuration plus complexe avec des données de test réelles");
    }
}

TEST_CASE("Test de follow_edges", "[follow_edges]") {
    // Note: Ces tests nécessitent une configuration complexe avec PackedGraph
    // et SnarlDistanceIndex. Il faudrait créer des mocks ou utiliser des
    // données de test réelles.
    SECTION("Documentation du test") {
        WARN("Les tests de follow_edges nécessitent une configuration plus complexe avec des données de test réelles");
    }
}

TEST_CASE("Test de fill_pretty_paths", "[fill_pretty_paths]") {
    // Note: Ces tests nécessitent une configuration complexe.
    // Il faudrait créer des mocks ou utiliser des données de test réelles.
    SECTION("Documentation du test") {
        WARN("Les tests de fill_pretty_paths nécessitent une configuration plus complexe avec des données de test réelles");
    }
}