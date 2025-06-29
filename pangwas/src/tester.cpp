#include "tester.hpp"
#include "utils.hpp"
#include <cmath>
#include <cassert>
#include <random>

//#define DEBUG_TESTER

using namespace std;
namespace pangwas {

FishersTester::FishersTester(const std::set<std::string>& samples_of_interest, double p_value_threshold, size_t sample_count) : 
    Tester(samples_of_interest), 
    p_value_threshold(p_value_threshold),
    total_sample_count(sample_count),
    cached_factorials(std::vector<double>(sample_count+1, std::numeric_limits<double>::max())) {}

bool FishersTester::is_associated(const std::set<std::string>& samples) {
    // For a contingency table a b
    //                         c d
    // p = (a+b)!(c+d)!(a+c)!(b+d)! / (a!b!c!d!n!)
    // 
    // a = # associated and on this allele
    // b = # unassociated and on this allele
    // c = # associated and not on this allele
    // d = # unassociated and not on this allele
    size_t a = 0;
    for (const std::string& sample : samples) {
        if (samples_of_interest.count(sample) != 0) {
            a += 1;
        }
    }
    size_t b = samples.size() - a;

    size_t c = samples_of_interest.size() - a;

    size_t d = total_sample_count - samples_of_interest.size() - b;

    return fishers_p_value(a, b, c, d) <= p_value_threshold;
}

double FishersTester::fishers_p_value(size_t a, size_t b, size_t c, size_t d) {
    #ifdef DEBUG_TESTER
   std::cerr << "Get fishers p value for values " << a << " " << b << " " << c << " " << d << std::endl;
    #endif
 
    double prob = log_fishers_probability(a, b, c, d);   

    double p = 0;

    // Since we only want the two ends, go up and down and skip the middle

    //The first tested a value that exceeded the probability going up
    size_t max_tested_a = 0;

    //The range of values of a that we can test given the marginals
    size_t max_a = std::min(a+b, a+c);
    size_t min_a = a >= d ? a-d : 0;
    #ifdef DEBUG_TESTER
   std::cerr << "\ttest values between " << min_a << " " << max_a << std::endl;
    #endif

    for (size_t test_a = min_a ; test_a <= max_a ; test_a++) {
        #ifdef DEBUG_TESTER
        assert(test_a >= 0);
        assert(a+b-test_a >= 0);
        assert(a+c-test_a >= 0);
        assert(d-a+test_a >= 0);
        assert(test_a + a + b - test_a + a + c - test_a + d - a + test_a == a+b+c+d);
        #endif

        double new_p = log_fishers_probability(test_a, a+b-test_a, a+c-test_a, d-a+test_a); 
        max_tested_a = test_a;

        // If this probability is less than or equal to (with some wiggle room) the probability of the original table
        // TODO: idk how much wiggle room
        if (new_p <= prob || is_equal(new_p, prob, std::numeric_limits<double>::epsilon()*20)) { 
            p += std::exp(new_p);
        } else {
            break;
        }
    }
    //Now go down
    for (int test_a = max_a ; test_a > max_tested_a ; test_a--) {
        #ifdef DEBUG_TESTER
        assert(test_a >= 0);
        assert(a+b-test_a >= 0);
        assert(a+c-test_a >= 0);
        assert(d-a+test_a >= 0);
        assert(test_a + a + b - test_a + a + c - test_a + d - a + test_a == a+b+c+d);
        #endif

        double new_p = log_fishers_probability(test_a, a+b-test_a, a+c-test_a, d-a+test_a); 
        if (new_p <= prob || is_equal(new_p, prob, std::numeric_limits<double>::epsilon()*20)) { 
            p += std::exp(new_p);
        } else {
            break;
        }
    }

    #ifdef DEBUG_TESTER
   std::cerr << "fishers p value: " << p << std::endl;
    #endif
    return p;
}

double FishersTester::log_fishers_probability(size_t a, size_t b, size_t c, size_t d) {
    #ifdef DEBUG_TESTER
   std::cerr << "\tfishers probability: " << a << " " << b << " " << c << " " << d << std::endl;
    #endif

    if (cached_probability.count(std::make_pair(std::make_pair(a, b), std::make_pair(c, d))) != 0) {
        #ifdef DEBUG_TESTER
       std::cerr << "\t\t" <<  cached_probability[std::make_pair(std::make_pair(a, b), std::make_pair(c, d))] << std::endl;
        #endif
        return cached_probability[std::make_pair(std::make_pair(a, b), std::make_pair(c, d))];
    }

    double p = log_factorial(a + b) + log_factorial(c + d) + log_factorial(a + c) + log_factorial(b + d)
              - log_factorial(a) - log_factorial(b) - log_factorial(c) - log_factorial(d) - log_factorial(a+b+c+d);

    cached_probability[std::make_pair(std::make_pair(a, b), std::make_pair(c, d))] = p;
    #ifdef DEBUG_TESTER
   std::cerr << "\t\t" <<  p << std::endl;
    #endif

    return p;
}

double FishersTester::log_factorial(size_t x) {
    //This shouldn't really happen but just in case, and for tests
    if ( x >= cached_factorials.size()) {
        cached_factorials.resize(x+1, std::numeric_limits<double>::max());
    }
    if (cached_factorials[x] != std::numeric_limits<double>::max()) {
        return cached_factorials[x];
    }
    double f = 0;
    for (size_t i = 1 ; i <= x ; i++) {
        f += std::log(i);
    }
    cached_factorials[x] = f;
    return f;
}

Chi2Tester::Chi2Tester(const std::set<std::string>& samples_of_interest, double p_value_threshold, size_t sample_count) : 
    Tester(samples_of_interest), 
    p_value_threshold(p_value_threshold),
    total_sample_count(sample_count),
    chi_squared_dist(1) {}

bool Chi2Tester::is_associated(const std::set<std::string>& samples) {
    // For a contingency table a b
    //                         c d
    //
    // a = # associated and on this allele
    // b = # unassociated and on this allele
    // c = # associated and not on this allele
    // d = # unassociated and not on this allele
    size_t a = 0;
    for (const std::string& sample : samples) {
        if (samples_of_interest.count(sample) != 0) {
            a += 1;
        }
    }

    size_t b = samples.size() - a;

    size_t c = samples_of_interest.size() - a;

    size_t d = total_sample_count - samples_of_interest.size() - b;

    // Get the p-value from chi-squared distribution
    #ifdef DEBUG_TESTER
    if (a != 0) {
        double p = p_value(a, b, c, d);
       std::cerr << "For counts: " << a << "\t" << b << endl
             << "            " << c << "\t" << d << std::endl;
       std::cerr << "\tChi2 p-value: " << p << std::endl; 
    }
    #endif
         
    return p_value(a, b, c, d) < p_value_threshold;
}

double Chi2Tester::p_value(size_t a, size_t b, size_t c, size_t d) {
    // TODO: It's probably better to use a lookup table but I can't find one
    return 1-boost::math::cdf(chi_squared_dist, test_statistic(a, b, c, d));
}


double Chi2Tester::test_statistic(size_t a, size_t b, size_t c, size_t d) {
    // For a contingency table a b
    //                         c d
    size_t total = a+b+c+d;

    if (total == 0) {
        return std::numeric_limits<double>::max();
    }

    double expected_a = (double)(a+b) * (a+c) / total; 
    double expected_b = (double)(a+b) * (b+d) / total;
    double expected_c = (double)(a+c) * (c+d) / total;
    double expected_d = (double)(b+d) * (c+d) / total;  

    if (expected_a == 0 || expected_b == 0  || expected_d == 0 || expected_d == 0) {
        return std::numeric_limits<double>::max();
    }

    double test_statistic = 0;
    test_statistic += std::pow((double)a - expected_a, 2) / expected_a;
    test_statistic += std::pow((double)b - expected_b, 2) / expected_b;
    test_statistic += std::pow((double)c - expected_c, 2) / expected_c;
    test_statistic += std::pow((double)d - expected_d, 2) / expected_d;


    return test_statistic;
}
   
}// end namespace pangwas
