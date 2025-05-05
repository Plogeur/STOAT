#include <catch2/catch_test_macros.hpp>
#include "../../src/list_snarl_paths.hpp"

using namespace std;
using namespace bdsg;
using handlegraph::step_handle_t;
using handlegraph::handle_t;
using handlegraph::net_handle_t;

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
        auto types = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "T");  // SNP
        REQUIRE(types[1] == "2");  // INS
        REQUIRE(types[2] == "0");  // DEL
    }

    SECTION("Cas SNP uniquement") {
        std::vector<std::tuple<std::string, size_t, size_t, size_t>> list_paths = {
            {"T", 0, 0, 3},
            {"G", 0, 0, 3}
        };
        auto types = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 2);
        REQUIRE(types[0] == "T");
        REQUIRE(types[1] == "G");
    }

    SECTION("Cas complexe") {
        std::vector<std::tuple<std::string, size_t, size_t, size_t>> list_paths = {
            {"_", 10, 20, 10},   // Complex path > 3
            {"TTTT", 0, 0, 3},   // Insertion
            {"", 0, 0, 2}        // Deletion
        };
        auto types = calcul_pos_type_variant(list_paths);
        REQUIRE(types.size() == 3);
        REQUIRE(types[0] == "CPX:10/20");
        REQUIRE(types[1] == "4");
        REQUIRE(types[2] == "0");
    }
}

TEST_CASE("Test simulated case", "[Path]") {

    std::unordered_map<std::string, std::vector<std::tuple<string, vector<string>, size_t, size_t, vector<string>>>> snarls_chr;
    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<bdsg::PackedGraph> pg;
    handlegraph::net_handle_t root;
    std::unique_ptr<bdsg::PackedPositionOverlay> pp_overlay;

    size_t children_threshold = 50;
    size_t path_length_threshold = 10000;
    std::unordered_set<std::string> ref_chr = {"ref"};
    bool only_snarl_parsing = false;
    string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
    string output_file = output_dir + "/snarl_analyse.tsv";

    SECTION("Simple SNP") {
        std::string pg_path = "../tests/graph_test/simple_snp.pg";
        std::string dist_path = "../tests/graph_test/simple_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "5_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5",">2>4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"C", "G"});
    }

    SECTION("3th SNP") {
        std::string pg_path = "../tests/graph_test/3th_snp.pg";
        std::string dist_path = "../tests/graph_test/3th_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "6_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>6",">2>4>6",">2>5>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"C", "T", "G"});
    }

    SECTION("4th") {
        std::string pg_path = "../tests/graph_test/4th.pg";
        std::string dist_path = "../tests/graph_test/4th.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "2_7");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5>7", ">2>4>6>7" ,">2>3>6>7"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 13);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"CPX:4/4", "CPX:5/5", "CPX:6/6"});
    }

    SECTION("deletion SNP") {
        std::string pg_path = "../tests/graph_test/deletion_snp.pg";
        std::string dist_path = "../tests/graph_test/deletion_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "4_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>4", ">2>3>4"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0", "C"});
    }

    SECTION("insert_deletion") {
        std::string pg_path = "../tests/graph_test/insert_deletion.pg";
        std::string dist_path = "../tests/graph_test/insert_deletion.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "4_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>4", ">2>3>4"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 12);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0", "3"});
    }

    SECTION("insert_snp") {
        std::string pg_path = "../tests/graph_test/insert_snp.pg";
        std::string dist_path = "../tests/graph_test/insert_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "5_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5", ">2>4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"C", "3"});
    }

    SECTION("inversion") {
        std::string pg_path = "../tests/graph_test/inversion.pg";
        std::string dist_path = "../tests/graph_test/inversion.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        
        // {chr : matrix(snarl, paths, start_pos, end_pos, type)}
        snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);        
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "4_1");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">1>2>3>4"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 12);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"C", "3"});
    }
}

std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, only_snarl_parsing);

