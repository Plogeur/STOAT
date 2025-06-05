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

TEST_CASE("Test simulated case", "[Path]") {

    std::unique_ptr<bdsg::SnarlDistanceIndex> stree;
    std::unique_ptr<bdsg::PackedGraph> pg;
    handlegraph::net_handle_t root;
    std::unique_ptr<bdsg::PackedPositionOverlay> pp_overlay;

    size_t children_threshold = 50;
    size_t path_length_threshold = 10000;
    std::unordered_set<std::string> ref_chr = {"ref"};
    bool only_snarl_parsing = false;
    std::string output_dir = "../tests/graph_test";
    string output_snarl_not_analyse = output_dir + "/snarl_not_analyse.tsv";
    string output_file = output_dir + "/snarl_analyse.tsv";

    SECTION("simple_snp") {
        std::string pg_path = "../tests/graph_test/simple_snp.pg";
        std::string dist_path = "../tests/graph_test/simple_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "5_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5",">2>4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"1", "1"});
    }

    SECTION("3th SNP") {
        std::string pg_path = "../tests/graph_test/3th_snp.pg";
        std::string dist_path = "../tests/graph_test/3th_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "6_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>6",">2>4>6",">2>5>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"1", "1", "1"});
    }

    SECTION("4th") {
        std::string pg_path = "../tests/graph_test/4th.pg";
        std::string dist_path = "../tests/graph_test/4th.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "2_7");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5>7", ">2>4>6>7" ,">2>3>6>7"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 13);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"4", "6", "5"});
    }

    SECTION("deletion_snp") {
        std::string pg_path = "../tests/graph_test/deletion_snp.pg";
        std::string dist_path = "../tests/graph_test/deletion_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "4_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>4", ">2>3>4"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0", "1"});
    }

    SECTION("insert_deletion") {
        std::string pg_path = "../tests/graph_test/insert_deletion.pg";
        std::string dist_path = "../tests/graph_test/insert_deletion.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
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
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "5_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5", ">2>4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"1", "3"});
    }

    SECTION("inversion") {
        std::string pg_path = "../tests/graph_test/inversion.pg";
        std::string dist_path = "../tests/graph_test/inversion.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "2_6");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>6", ">2>3>*>5>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 15);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0", "6/6"});
        
        REQUIRE(std::get<1>(snarls_chr["ref"][1]) == std::vector<std::string>{">3>4>5", ">3<4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][1]) == 9);
        REQUIRE(std::get<3>(snarls_chr["ref"][1]) == 12);
        REQUIRE(std::get<4>(snarls_chr["ref"][1]) == std::vector<std::string>{"2", "2"});
    }

    SECTION("large_del") {
        std::string pg_path = "../tests/graph_test/large_del.pg";
        std::string dist_path = "../tests/graph_test/large_del.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        
        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 3);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "2_9");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>9",">2>3>*>8>9"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 9);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0", "9/10"}); // correct CPX by boundary the first and end complexe chain node not the snarl boundary

        REQUIRE(std::get<0>(snarls_chr["ref"][1]) == "6_8");
        REQUIRE(std::get<1>(snarls_chr["ref"][1]) == std::vector<std::string>{">6>8",">6>7>8"});
        REQUIRE(std::get<2>(snarls_chr["ref"][1]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][1]) == 9);
        REQUIRE(std::get<4>(snarls_chr["ref"][1]) == std::vector<std::string>{"0", "1"});

        REQUIRE(std::get<0>(snarls_chr["ref"][2]) == "3_6");
        REQUIRE(std::get<1>(snarls_chr["ref"][2]) == std::vector<std::string>{">3>5>6",">3>4>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][2]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][2]) == 9);
        REQUIRE(std::get<4>(snarls_chr["ref"][2]) == std::vector<std::string>{"1", "1"});
    }

    SECTION("linear") {
        std::string pg_path = "../tests/graph_test/linear.pg";
        std::string dist_path = "../tests/graph_test/linear.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);
        REQUIRE(snarls_chr.size() == 0);
    }

    SECTION("loop_simple") {
        std::string pg_path = "../tests/graph_test/loop_simple.pg";
        std::string dist_path = "../tests/graph_test/loop_simple.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 1, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "5_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>5", ">2>3>3>5", ">2>4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"1", "2", "2"});
    }

    SECTION("loop") {
        std::string pg_path = "../tests/graph_test/loop.pg";
        std::string dist_path = "../tests/graph_test/loop.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 2, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "5_1");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">1>2>3>5", ">1>2>3>2>3>5", ">1>2>3>2>3>2>3>5", ">1>2>3>2>3>2>4>5", ">1>2>3>2>4>5", ">1>2>4>5"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 4);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"5", "10", "15", "16", "11", "6"});
    }

    SECTION("loop_double") {
        std::string pg_path = "../tests/graph_test/loop_double.pg";
        std::string dist_path = "../tests/graph_test/loop_double.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 2, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "7_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>4>5>7", ">2>3>4>5>3>4>5>7", ">2>3>4>5>3>4>5>3>4>5>7", ">2>3>4>5>3>4>3>4>5>7", ">2>3>4>3>4>5>7", ">2>3>4>3>4>5>3>4>5>7", ">2>3>4>3>4>3>4>5>7", ">2>6>7"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 11);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"6", "12", "18", "16", "10", "16", "14", "2"});
    }

    SECTION("loop_plus") {
        std::string pg_path = "../tests/graph_test/loop_plus.pg";
        std::string dist_path = "../tests/graph_test/loop_plus.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 1, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "8_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>*>6>8",">2>3>*>6>3>*>6>8",">2>7>8"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"3/8","3/8","1"});

        REQUIRE(snarls_chr["ref"].size() == 2);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "3_6");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>*>6>8",">2>3>*>6>3>*>6>8",">2>7>8"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"3/8","3/8","1"});

    }

    SECTION("repetition") {
        std::string pg_path = "../tests/graph_test/repetition.pg";
        std::string dist_path = "../tests/graph_test/repetition.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "6_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>6", ">2>3>6", ">2>3>4>6", ">2>3>4>5>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 9);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0", "3", "6", "9"});
    }

    SECTION("complex_ins") {
        std::string pg_path = "../tests/graph_test/complex_ins.pg";
        std::string dist_path = "../tests/graph_test/complex_ins.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 1);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "8_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>8", ">2>3>4>6>8",">2>3>5>6>8",">2>3>5>7>8",">2>7>8"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0","3","3","3","1"});
    }

    SECTION("snp_and_nested_snp") {
        std::string pg_path = "../tests/graph_test/snp_and_nested_snp.pg";
        std::string dist_path = "../tests/graph_test/snp_and_nested_snp.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "8_2");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>3>*>6>8",">2>7>8"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"3/4","1"});

        REQUIRE(std::get<0>(snarls_chr["ref"][1]) == "6_3");
        REQUIRE(std::get<1>(snarls_chr["ref"][1]) == std::vector<std::string>{">3>4>6",">3>5>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][1]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][1]) == 10);
        REQUIRE(std::get<4>(snarls_chr["ref"][1]) == std::vector<std::string>{"1","2"});
    }

    SECTION("nested_plus") {
        std::string pg_path = "../tests/graph_test/nested_plus.pg";
        std::string dist_path = "../tests/graph_test/nested_plus.dist";

        std::tie(stree, pg, root, pp_overlay) = parse_graph_tree(pg_path, dist_path);
        auto snarls = save_snarls(*stree, root, *pg, ref_chr, *pp_overlay);
        auto snarls_chr = loop_over_snarls_write(*stree, snarls, *pg, output_file, output_snarl_not_analyse, children_threshold, path_length_threshold, 0, only_snarl_parsing);

        REQUIRE(snarls_chr.size() == 1);
        REQUIRE(snarls_chr["ref"].size() == 2);
        REQUIRE(std::get<0>(snarls_chr["ref"][0]) == "2_8");
        REQUIRE(std::get<1>(snarls_chr["ref"][0]) == std::vector<std::string>{">2>8", ">2>3>*>6>7>8", ">2>3>*>6>8"});
        REQUIRE(std::get<2>(snarls_chr["ref"][0]) == 8);
        REQUIRE(std::get<3>(snarls_chr["ref"][0]) == 13);
        REQUIRE(std::get<4>(snarls_chr["ref"][0]) == std::vector<std::string>{"0","5/5", "4/4"});

        REQUIRE(std::get<0>(snarls_chr["ref"][1]) == "3_6");
        REQUIRE(std::get<1>(snarls_chr["ref"][1]) == std::vector<std::string>{">3>5>6", ">3>4>6"});
        REQUIRE(std::get<2>(snarls_chr["ref"][1]) == 9);
        REQUIRE(std::get<3>(snarls_chr["ref"][1]) == 12);
        REQUIRE(std::get<4>(snarls_chr["ref"][1]) == std::vector<std::string>{"2","2"});
    }
}
