#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../path_association_finder.hpp"


namespace pangwas{

class TestPathAssociationFinder : PathAssociationFinder {
    public: 
    TestPathAssociationFinder(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                           std::string test, const std::set<std::string>& samples_of_interest, std::string reference_name,
                           string output_format, std::ostream& out_associated, std::ostream& out_unassociated,
                           size_t allele_size_limit, double p_value) :
        PathAssociationFinder(graph, distance_index, test, samples_of_interest, reference_name, output_format, out_associated, 
                          out_unassociated, allele_size_limit, p_value) {} 
    using PathAssociationFinder::get_walk_sets;
    using PathAssociationFinder::get_start_edge_sets;
    using PathAssociationFinder::is_snarl_associated;
};

TEST_CASE( "Path association finder one node",
          "[path_finder]" ) {


    bdsg::HashGraph graph;
        
    handlegraph::handle_t n1 = graph.create_handle("GCAAACAGATT");

    handlegraph::path_handle_t path = graph.create_path_handle("path");
    graph.append_step(path, n1);

    // vg isn't included so the distance index can only be built from the command line
    graph.serialize("test.hg");
    int built = system("vg index -j test.dist test.hg"); 
    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("test.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);

    SECTION("Make association finder") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing
        TestPathAssociationFinder af(*path_graph, distance_index, "exact", std::set<std::string>(), "a", "b", std::cout, std::cout, 10, 0.0);
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}
TEST_CASE( "Path association finder nested bubbles",
          "[path_finder]" ) {

    /*
                       5
                     /   \
            1       4 ----6    8
          /   \   /         \ / \
        0       3  ----------7---9
          \   /
            2

   */

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = { "C", "C", "C", "A", "T", "C", "A", "C", "A", "A"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(graph.create_handle(seq));
    }

    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[0], nodes[2]);
    graph.create_edge(nodes[1], nodes[3]);
    graph.create_edge(nodes[2], nodes[3]);
    graph.create_edge(nodes[3], nodes[4]);
    graph.create_edge(nodes[3], nodes[7]);
    graph.create_edge(nodes[4], nodes[5]);
    graph.create_edge(nodes[4], nodes[6]);
    graph.create_edge(nodes[5], nodes[6]);
    graph.create_edge(nodes[6], nodes[7]);
    graph.create_edge(nodes[7], nodes[8]);
    graph.create_edge(nodes[7], nodes[9]);
    graph.create_edge(nodes[8], nodes[9]);

    // TODO one of these should really be the reference but idk how to add reference paths to a graph
    std::vector<std::vector<std::size_t>> paths_seqs = { {0, 1, 3, 4, 5, 6, 7}, {0, 1, 3, 4, 6, 7}, {0, 2, 3, 7}, {0, 2, 3, 4, 6, 7}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < paths_seqs.size() ; path_i++) {
        paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
        for (size_t node_i : paths_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    // vg isn't included so the distance index can only be built from the command line
    graph.serialize("test.hg");
    int built = system("vg index -j test.dist test.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("test.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[1], &graph)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[4], &graph)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[5], &graph)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[8], &graph)));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);
    handlegraph::net_handle_t nested_chain = distance_index.get_parent(snarl3);

    // snarl3 should be associated
    std::set<std::string> samples ({"path1", "path3"});
    TestPathAssociationFinder af(*path_graph, distance_index, "exact", samples, "a", "b", std::cout, std::cout, 10, 0.0);


    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} and {2,3}
        std::vector<std::set<pangwas::sample_hap_t>> walks1 = af.get_walk_sets(snarl1);
        REQUIRE(walks1.size() == 2);
        for ( const auto& set : walks1) {
            REQUIRE(set.size() == 2);
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2]), pangwas::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }

        // Should be {0,1,3} and {2}
        std::vector<std::set<pangwas::sample_hap_t>> walks2 = af.get_walk_sets(snarl2);
        REQUIRE(walks2.size() == 2);
        for ( const auto& set : walks2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[1]), pangwas::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }

        // Should be {0}, {1,3} and {2}
        std::vector<std::set<pangwas::sample_hap_t>> walks3 = af.get_walk_sets(snarl3);
        REQUIRE(walks3.size() == 3);
        for ( const auto& set : walks3) {
            REQUIRE(((set.size() == 2) || (set.size() == 1)));
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0])}) ) ||
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[1]), pangwas::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }
    }
    SECTION("get start edge sets") {
        // Should be {0,1} and {2,3}
        std::vector<std::set<pangwas::sample_hap_t>> edges1 = af.get_start_edge_sets(snarl1);
        REQUIRE(edges1.size() == 2);
        for ( const auto& set : edges1) {
            REQUIRE(set.size() == 2);
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2]), pangwas::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }

        // Should be {0,1,3} and {2}
        std::vector<std::set<pangwas::sample_hap_t>> edges2 = af.get_start_edge_sets(snarl2);
        REQUIRE(edges2.size() == 2);
        for ( const auto& set : edges2) {
            REQUIRE(((set.size() == 3) || (set.size() == 1)));
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[1]), pangwas::get_sample_and_haplotype(*path_graph, paths[3])})) || 
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }

        // Should be {0} and {1,3}
        std::vector<std::set<pangwas::sample_hap_t>> edges3 = af.get_start_edge_sets(snarl3);
        REQUIRE(edges3.size() == 2);
        for ( const auto& set : edges3) {
            REQUIRE(((set.size() == 2) || (set.size() == 1)));
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0])}) ) ||
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[1]), pangwas::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("test snarl_is_associated") {
        REQUIRE(!af.is_snarl_associated(snarl1).first);
        REQUIRE(!af.is_snarl_associated(snarl2).first);
        REQUIRE(af.is_snarl_associated(snarl3).first);
        REQUIRE(af.is_snarl_associated(snarl3).second.size() == 3);
        REQUIRE(!af.is_snarl_associated(snarl4).first);
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}

TEST_CASE( "Path association finder looping snarl",
          "[path_finder]" ) {

    /*

             --------
            |   2    |
            \ / \    /
        0 ---1---3--4----5

    */

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C", "T",  "AAAAAAAAA"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(graph.create_handle(seq));
    }

    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[1], nodes[2]);
    graph.create_edge(nodes[1], nodes[3]);
    graph.create_edge(nodes[2], nodes[3]);
    graph.create_edge(nodes[3], nodes[4]);
    graph.create_edge(nodes[4], nodes[1]);
    graph.create_edge(nodes[4], nodes[5]);


    // Paths 0 and 2 take the insertion, but paths 1 and 2 take the duplication, and the deletion
    std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 2, 3, 4, 5}, {0, 1, 3, 4, 1, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
        paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
        for (size_t node_i : path_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    // vg isn't included so the distance index can only be built from the command line
    graph.serialize("test.hg");
    int built = system("vg index -j test.dist test.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("test.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[2], &graph)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    std::set<std::string> samples ({"path1", "path2"});
    TestPathAssociationFinder af(*path_graph, distance_index, "exact", samples, "a", "b", std::cout, std::cout, 10, 0.0);

    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0} and {1,2}
        std::vector<std::set<pangwas::sample_hap_t>> walks1 = af.get_walk_sets(snarl1);
        REQUIRE(walks1.size() == 2);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[1]), pangwas::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0])}))));
        }

        // Should be {0,2} and {1}
        std::vector<std::set<pangwas::sample_hap_t>> walks2 = af.get_walk_sets(snarl2);
        REQUIRE(walks2.size() == 2);
        for ( const auto& set : walks2) {
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[1])}))));
        }
    }
    SECTION("get_start_edge_set") {

        // Should be {0,22 and {1,2}
        std::vector<std::set<pangwas::sample_hap_t>> edges2 = af.get_start_edge_sets(snarl2);
        REQUIRE(edges2.size() == 2);
        for ( const auto& set : edges2) {
            REQUIRE(set.size() == 2);
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                      (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[1]), pangwas::get_sample_and_haplotype(*path_graph, paths[2])}))));
        }
    }
    SECTION("test snarl_is_associated") {
        // idk about this but the duplication is flagged and also the deletion
        REQUIRE(af.is_snarl_associated(snarl1).first);
        REQUIRE(af.is_snarl_associated(snarl1).second.size() == 3);
        REQUIRE(af.is_snarl_associated(snarl2).first);
        REQUIRE(af.is_snarl_associated(snarl2).second.size() == 3);
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}
TEST_CASE( "Path association finder bubble with three nodes",
          "[path_finder]" ) {

    /*
           1    
         /   \
        0--2--4
         \   /
           3

    */

    bdsg::HashGraph graph;

    std::vector<std::string> sequences = {"AAAAAAAAAA", "A", "G", "C",  "AAAAAAAAA"};

    std::vector<handlegraph::handle_t> nodes;
    for (auto& seq : sequences) {
        nodes.emplace_back(graph.create_handle(seq));
    }

    graph.create_edge(nodes[0], nodes[1]);
    graph.create_edge(nodes[0], nodes[2]);
    graph.create_edge(nodes[0], nodes[3]);
    graph.create_edge(nodes[1], nodes[4]);
    graph.create_edge(nodes[2], nodes[4]);
    graph.create_edge(nodes[3], nodes[4]);


    // Two paths go through node 2, path 2 is associated
    std::vector<std::vector<std::size_t>> path_seqs = { {0, 1, 4}, {0, 1, 4}, {0, 2, 4}, {0, 3, 4}};
    std::vector<handlegraph::path_handle_t> paths;

    for (int path_i = 0 ; path_i < path_seqs.size() ; path_i++) {
        paths.emplace_back(graph.create_path_handle("path"+std::to_string(path_i)));
        for (size_t node_i : path_seqs[path_i]) {
            graph.append_step(paths.back(), nodes[node_i]);
        }
    }

    // vg isn't included so the distance index can only be built from the command line
    graph.serialize("test.hg");
    int built = system("vg index -j test.dist test.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("test.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;
    auto path_graph = overlay_helper.apply(&graph);


    handlegraph::net_handle_t snarl = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[2], &graph)));


    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    std::set<std::string> samples ({"path2"});
    TestPathAssociationFinder af(*path_graph, distance_index, "exact", samples, "a", "b", std::cout, std::cout, 10, 0.0);

    SECTION("get_walk_set") {
        // This isn't really a good test because all the snarls are regular

        // Should be {0,1} {2} {3}
        std::vector<std::set<pangwas::sample_hap_t>> walks1 = af.get_walk_sets(snarl);
        REQUIRE(walks1.size() == 3);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("get_start_edge_set") {

        // Should be {0,1} {2} {3}
        std::vector<std::set<pangwas::sample_hap_t>> walks1 = af.get_walk_sets(snarl);
        REQUIRE(walks1.size() == 3);
        for ( const auto& set : walks1) {
            REQUIRE( ((set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[0]), pangwas::get_sample_and_haplotype(*path_graph, paths[1])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[2])})) || 
                     (set == std::set<pangwas::sample_hap_t> ({pangwas::get_sample_and_haplotype(*path_graph, paths[3])}))));
        }
    }
    SECTION("test snarl_is_associated") {
        REQUIRE(af.is_snarl_associated(snarl).first);
        // This should have one sample from each set
        REQUIRE(af.is_snarl_associated(snarl).second.size() == 3);
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}
}
