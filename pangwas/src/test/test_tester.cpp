#include <catch.hpp>
#include <iostream>
#include <iomanip>
#include "../utils.hpp"
#include "../tester.hpp"

namespace pangwas {

class TestFishersTester : FishersTester {
    public:
    TestFishersTester(const std::set<std::string>& samples_of_interest, double p_cutoff, size_t sample_count) :
        FishersTester(samples_of_interest, p_cutoff, sample_count) {}
    using FishersTester::log_fishers_probability;
    using FishersTester::fishers_p_value;
    using FishersTester::is_associated;
};
class TestChi2Tester : Chi2Tester {
    public:
    TestChi2Tester(const std::set<std::string>& samples_of_interest, double p_cutoff, size_t sample_count) :
        Chi2Tester(samples_of_interest, p_cutoff, sample_count) {}
    using Chi2Tester::test_statistic;
    using Chi2Tester::p_value;
};

TEST_CASE( "Exact test ","[tester]" ) {
    std::set<std::string> samples_of_interest ({"1", "2", "3"});
    ExactTester tester(samples_of_interest);


    SECTION("Test tester") {
        REQUIRE(tester.is_associated(std::set<std::string>({"2", "1", "3"})));
        REQUIRE(!tester.is_associated(std::set<std::string>({"2", "1"})));
        REQUIRE(!tester.is_associated(std::set<std::string>({"2", "1", "3", "4"})));
        REQUIRE(!tester.is_associated(std::set<std::string>({"2", "1", "4"})));
    }
}
TEST_CASE( "Fishers test ","[tester][fishers]" ) {
    std::set<std::string> samples_of_interest ({"1", "2", "3", "4"});
    TestFishersTester tester(samples_of_interest, 0.05, 30);


    SECTION("Test fishers probability") {
        REQUIRE(tester.is_associated(std::set<std::string>({"2", "1", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"})));
        REQUIRE(tester.is_associated(std::set<std::string>({"2", "1", "3", "4"})));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(1, 9, 11, 3)),
                 0.001346076, 0.000001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(0, 10, 12, 2)),
                 0.000033652, 0.000001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(4, 10, 0, 16)),
                 0.036, 0.001));

        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(2, 6, 5, 0)),
                 0.0163, 0.0001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(3, 5, 4, 1)),
                 0.163, 0.001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(4, 4, 3, 2)),
                 0.408, 0.001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(5, 3, 2, 3)),
                 0.326, 0.001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(6, 2, 1, 4)),
                 0.0816, 0.0001));
        REQUIRE(is_equal(std::exp(tester.log_fishers_probability(7, 1, 0, 5)),
                 0.00466, 0.0001));
    }

    SECTION("Test fishers p-value") {
        REQUIRE(is_equal(tester.fishers_p_value(1, 9, 11, 3),
                 0.00275946, 0.000001));
        REQUIRE(is_equal(tester.fishers_p_value(6, 2, 1, 4),
                 0.10256410256410257, 0.00000001));
        REQUIRE(is_equal(tester.fishers_p_value(6, 1, 2, 7),
                 0.0406, 0.0001));
        REQUIRE(is_equal(tester.fishers_p_value(0, 4, 4, 7),
                 0.5165, 0.0001));
        REQUIRE(is_equal(tester.fishers_p_value(5, 5, 5, 5),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(5, 0, 8, 3),
                 0.5089, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(5, 0, 0, 5),
                 0.0079, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(10, 0, 0, 10),
                 0.0001, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(10, 10, 0, 0),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(10, 0, 10, 0),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(0, 10, 10, 0),
                 0.0001, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(0, 10, 0, 10),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(1, 9, 0, 10),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(9, 1, 1, 9),
                 0.0011, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(7, 2, 2, 9),
                 0.0216, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(7, 6, 4, 3),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.fishers_p_value(7, 1, 4, 8),
                 0.0281, 0.001));
    }

}
TEST_CASE( "Chi2 test ","[tester][chi2]" ) {
    std::set<std::string> samples_of_interest ({"1", "2", "3", "4"});
    TestChi2Tester tester(samples_of_interest, 0.01, 30);

    SECTION("Test chi2 test statistic") {
        REQUIRE(is_equal(tester.test_statistic(36, 14, 
                                               30, 25),
                 3.4177, 0.001));
        REQUIRE(is_equal(tester.test_statistic(50, 50, 
                                               50, 50),
                 0.0, 0.001));
        REQUIRE(is_equal(tester.test_statistic(50, 0, 
                                               50, 50),
                 37.5, 0.1));
        REQUIRE(is_equal(tester.test_statistic(50, 50, 
                                               0, 50),
                 37.5, 0.1));
        REQUIRE(is_equal(tester.test_statistic(0, 50, 
                                               50, 50),
                 37.5, 0.1));
        REQUIRE(is_equal(tester.test_statistic(50, 50, 
                                               50, 0),
                 37.5, 0.1));
        REQUIRE(is_equal(tester.test_statistic(0, 50, 
                                               50, 0),
                 100.0, 0.1));
        REQUIRE(is_equal(tester.test_statistic(50, 0, 
                                               0, 50),
                 100.0, 0.1));
        REQUIRE(is_equal(tester.test_statistic(345, 125, 
                                               234, 234),
                 54.3705, 0.001));
        REQUIRE(is_equal(tester.test_statistic(12, 65, 
                                               23, 76),
                 1.5901, 0.001));
        REQUIRE(is_equal(tester.test_statistic(65, 87, 
                                               23, 27),
                 0.1604, 0.001));
        REQUIRE(is_equal(tester.test_statistic(65, 23, 
                                               56, 82),
                 23.9312, 0.001));
        REQUIRE(is_equal(tester.test_statistic(0, 0, 
                                               0, 100),
                 std::numeric_limits<double>::max(), 0.001));
    }
    SECTION("Test chi2 p-value") {
        REQUIRE(is_equal(tester.p_value(36, 14, 
                                        30, 25),
                 0.0645, 0.001));
        REQUIRE(is_equal(tester.p_value(50, 50, 
                                        50, 50),
                 1.0, 0.001));
        REQUIRE(is_equal(tester.p_value(50, 0, 
                                        50, 50),
                 0.0, 0.1));
        REQUIRE(is_equal(tester.p_value(50, 50, 
                                        0, 50),
                 0.0, 0.1));
        REQUIRE(is_equal(tester.p_value(0, 50, 
                                        50, 50),
                 0.0, 0.1));
        REQUIRE(is_equal(tester.p_value(50, 50, 
                                        50, 0),
                 0.0, 0.1));
        REQUIRE(is_equal(tester.p_value(0, 50, 
                                        50, 0),
                 0.0, 0.1));
        REQUIRE(is_equal(tester.p_value(50, 0, 
                                        0, 50),
                 0.0, 0.1));
        REQUIRE(is_equal(tester.p_value(345, 125, 
                                        234, 234),
                 0.0, 0.001));
        REQUIRE(is_equal(tester.p_value(12, 65, 
                                        23, 76),
                 0.2073, 0.001));
        REQUIRE(is_equal(tester.p_value(65, 87, 
                                        23, 27),
                 0.6888, 0.001));
        REQUIRE(is_equal(tester.p_value(65, 23, 
                                        56, 82),
                 0.0, 0.001));
    }

}

} //end namespace pangwas
