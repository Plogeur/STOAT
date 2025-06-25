#include <catch.hpp>
#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../association_finder.hpp"
#include <regex>


namespace pangwas{

class TestAssociationFinder : AssociationFinder {
    public: 
    TestAssociationFinder(const handlegraph::PathPositionHandleGraph& graph, const bdsg::SnarlDistanceIndex& distance_index,
                            std::string test, const std::set<std::string>& samples_of_interest,std::string reference_name,
                           string output_format, std::ostream& out_associated, std::ostream& out_unassociated,
                           size_t allele_size_limit, double p_value) :
        AssociationFinder(graph, distance_index, test,samples_of_interest,  reference_name, output_format, out_associated, 
                          out_unassociated, allele_size_limit, p_value) {} 
    using AssociationFinder::path_range_t;
    using AssociationFinder::get_coordinates_of_snarl;
    using AssociationFinder::write_fasta_of_snarl;
    using AssociationFinder::graph;
    std::vector<std::set<std::string>> partition_samples_in_snarl(const handlegraph::net_handle_t& snarl) const {
        return std::vector<std::set<std::string>>();
    }
};

TEST_CASE( "Association finder one node",
          "[base_finder]" ) {


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

    SECTION("Make association finder") {
        // There isn't much to do with one node so just make sure we can run the constructor without crashing
        TestAssociationFinder af(*overlay_helper.apply(&graph), distance_index, "exact", std::set<std::string>(), "a", "b", std::cout, std::cout, 10, 0.0);
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}
TEST_CASE( "Association finder nested bubbles",
          "[base_finder]" ) {

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
    std::vector<std::vector<std::size_t>> paths = { {0, 1, 3, 4, 5, 6, 7}, {0, 1, 3, 4, 6, 7}, {0, 2, 3, 7}, {0, 2, 3, 4, 6, 7}};

    for (int path_i = 0 ; path_i < paths.size() ; path_i++) {
        handlegraph::path_handle_t path = graph.create_path_handle("path"+std::to_string(path_i));
        for (size_t node_i : paths[path_i]) {
            graph.append_step(path, nodes[node_i]);
        }
    }

    // vg isn't included so the distance index can only be built from the command line
    graph.serialize("test.hg");
    int built = system("vg index -j test.dist test.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("test.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;


    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[1], &graph)));
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[4], &graph)));
    handlegraph::net_handle_t snarl3 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[5], &graph)));
    handlegraph::net_handle_t snarl4 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[8], &graph)));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);
    handlegraph::net_handle_t nested_chain = distance_index.get_parent(snarl3);

    std::stringstream out;

    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    TestAssociationFinder af(*overlay_helper.apply(&graph), distance_index, "exact", std::set<std::string>(), "path0", "b", out, out, 10, 0.0);


    SECTION("Test get_coordinates_of_snarl for specific paths") {
        std::vector<TestAssociationFinder::path_range_t> ranges = af.get_coordinates_of_snarl(snarl1, false, "path0", false);
        REQUIRE(ranges.size() == 1);
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].start)) == "path0");
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].end)) == "path0");
        REQUIRE(af.graph.get_position_of_step(ranges[0].start) == 0);
        REQUIRE(af.graph.get_position_of_step(ranges[0].end) == 2);

        ranges = af.get_coordinates_of_snarl(snarl2, false, "path0", false);
        REQUIRE(ranges.size() == 1);
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].start)) == "path0");
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].end)) == "path0");
        REQUIRE(af.graph.get_position_of_step(ranges[0].start) == 2);
        REQUIRE(af.graph.get_position_of_step(ranges[0].end) == 6);

        ranges = af.get_coordinates_of_snarl(snarl3, false, "path0", false);
        REQUIRE(ranges.size() == 1);
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].start)) == "path0");
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].end)) == "path0");
        REQUIRE(af.graph.get_position_of_step(ranges[0].start) == 3);
        REQUIRE(af.graph.get_position_of_step(ranges[0].end) == 5);

        ranges = af.get_coordinates_of_snarl(snarl4, false, "path0", false);
        REQUIRE(ranges.size() == 0);
    }
    SECTION("Test get_coordinates_of_snarl for all paths") {
        std::vector<TestAssociationFinder::path_range_t> ranges = af.get_coordinates_of_snarl(snarl1, false, "", true);
        REQUIRE(ranges.size() == 4);
        for (const auto& range : ranges) {
            REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(range.start)) == 
                    graph.get_path_name(graph.get_path_handle_of_step(range.end)));
            REQUIRE(af.graph.get_position_of_step(range.start) == 0);
            REQUIRE(af.graph.get_position_of_step(range.end) == 2);
        }
    }
    SECTION("Test fasta output") {
        std::unordered_set<std::string> empty_set;
        af.write_fasta_of_snarl(snarl1, empty_set);
        std::string test;
        std::regex match(">snarl:1-4\\|path0:1-2\\|path[0-3]:1-2");
        while (std::getline(out, test)) {
            cerr << test << endl;
            if (test.substr(0,1) == ">") {
                REQUIRE(std::regex_match(test, match));
            } else {
                REQUIRE(test == "C");
            }
        }
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}

TEST_CASE( "Association finder looping snarl",
          "[base_finder]" ) {

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


    // Paths 0 and 2 take the insertion, but paths 1 and 2 take the duplication
    std::vector<std::vector<std::size_t>> paths = { {0, 1, 2, 3, 4, 5}, {0, 1, 3, 4, 1, 3, 4, 5}, {0, 1, 2, 3, 4, 1, 3, 4, 5}};

    for (int path_i = 0 ; path_i < paths.size() ; path_i++) {
        handlegraph::path_handle_t path = graph.create_path_handle("path"+std::to_string(path_i));
        for (size_t node_i : paths[path_i]) {
            graph.append_step(path, nodes[node_i]);
        }
    }

    // vg isn't included so the distance index can only be built from the command line
    graph.serialize("test.hg");
    int built = system("vg index -j test.dist test.hg"); 

    bdsg::SnarlDistanceIndex distance_index;
    distance_index.deserialize("test.dist");

    bdsg::PathPositionOverlayHelper overlay_helper;


    // Nested snarl
    handlegraph::net_handle_t snarl2 = distance_index.get_parent(distance_index.get_parent(distance_index.get_net(nodes[2], &graph)));
    // Duplication snarl
    handlegraph::net_handle_t snarl1 = distance_index.get_parent(distance_index.get_parent(snarl2));
    handlegraph::net_handle_t root_chain = distance_index.get_parent(snarl1);


    // This file is meant to test the base association finder but since it is technically an interface with some implementations,
    // build the path version and only test the base functions
    TestAssociationFinder af(*overlay_helper.apply(&graph), distance_index, "exact", std::set<std::string>(), "a", "b", std::cout, std::cout, 10, 0.0);


    SECTION("Test get_coordinates_of_snarl for specific paths") {
        std::vector<TestAssociationFinder::path_range_t> ranges = af.get_coordinates_of_snarl(snarl2, false, "path2", false);
        //Technically I think it isn't guaranteed to be in order
        REQUIRE(ranges.size() == 2);
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].start)) == "path2");
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[0].end)) == "path2");
        REQUIRE(af.graph.get_position_of_step(ranges[0].start) == 10);
        REQUIRE(af.graph.get_position_of_step(ranges[0].end) == 12);

        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[1].start)) == "path2");
        REQUIRE(graph.get_path_name(graph.get_path_handle_of_step(ranges[1].end)) == "path2");
        REQUIRE(af.graph.get_position_of_step(ranges[1].start) == 14);
        REQUIRE(af.graph.get_position_of_step(ranges[1].end) == 15);
    }
    SECTION("Test get_coordinates_of_snarl for all paths") {
        std::vector<TestAssociationFinder::path_range_t> ranges = af.get_coordinates_of_snarl(snarl2, false, "", true);
        REQUIRE(ranges.size() == 5);
    }

    // Remember to clean up the files made here
    int removed = system("rm -f test.hg test.dist"); 
}
}
