#ifndef PANGWAS_TESTER_HPP_INCLUDED
#define PANGWAS_TESTER_HPP_INCLUDED

#include <set>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;
namespace pangwas {
/***
    Template class for testing associations.
***/

class Tester {
    public:

    Tester (const std::set<std::string>& samples_of_interest) :
        samples_of_interest(samples_of_interest) 
        {}; 

    /// Is the set of samples associated? 
    virtual bool is_associated(const std::set<std::string>& samples) = 0;

    protected:

    // The samples known to be associated with the trait of interest
    // This should really be const but since I need a shared_ptr to it it can't be
    const std::set<std::string>& samples_of_interest;
 
};

/***
    Test that an set of samples exactly matches the set of samples of interest
***/
class ExactTester : public Tester { 
    public:

    ExactTester(const std::set<std::string>& samples_of_interest) :
        Tester(samples_of_interest) 
        {}

    bool is_associated(const std::set<std::string>& samples) {
        return samples == samples_of_interest;
    }
};

/***
    Fishers test
***/
class FishersTester : public Tester {
    public:

        FishersTester(const std::set<std::string>& samples_of_interest, double p_cutoff, size_t sample_count);

        bool is_associated(const std::set<std::string>& samples) ;

    protected:
        ///////////////////////// Helper functions

        // These are all really static functions, except that we cache the values
        // Get the log of the factorial of a number. Use cached factorials
        double log_factorial(size_t x);

        // The p value of a fishers test
        double fishers_p_value(size_t a, size_t b, size_t c, size_t d);

        // Get the log probability of the counts
        // a = # associated and on this allele
        // b = # unassociated and on this allele
        // c = # associated and not on this allele
        // d = # unassociated and not on this allele
        double log_fishers_probability(size_t a, size_t b, size_t c, size_t d);

        ///////////////////////// Member variables

        //What p value do we use to determine significance?
        double p_value_threshold = 0.05;

        // How many samples are there total
        size_t total_sample_count = 0;

        // Remember results that have already been calculated
        // Map pair of pairs <<a, b>, <c, d>> to the log probability
        std::map<std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t>>, double> cached_probability;
        std::vector<double> cached_factorials;
};

/***
    Chi-squared test
***/
class Chi2Tester : public Tester {
    public:

        Chi2Tester(const std::set<std::string>& samples_of_interest, double p_cutoff, size_t sample_count);

        bool is_associated(const std::set<std::string>& samples) ;

    protected:
        ///////////////////////// Helper functions

        // Get the chi-squared test statistic given counts
        // a = # associated and on this allele
        // b = # unassociated and on this allele
        // c = # associated and not on this allele
        // d = # unassociated and not on this allele
        double test_statistic(size_t a, size_t b, size_t c, size_t d);

        // Get the p-value. This is pulled out to test
        double p_value(size_t a, size_t b, size_t c, size_t d);

        ///////////////////////// Member variables

        //What p value do we use to determine significance?
        double p_value_threshold = 0.05;

        // How many samples are there total
        size_t total_sample_count = 0;

        // The chi-squared distribution for 1 degree of freedom
        boost::math::chi_squared chi_squared_dist;
};

}// end namespace pangwas
#endif
