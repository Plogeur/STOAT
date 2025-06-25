#include "utils.hpp"

using namespace std;
namespace pangwas {
std::string get_sample_name_from_path(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path) {

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        return graph.get_locus_name(path);
    } else {
        return graph.get_sample_name(path);
    }

}
sample_hap_t get_sample_and_haplotype(const handlegraph::PathPositionHandleGraph& graph, const handlegraph::path_handle_t& path) {
    sample_hap_t result;

    if (graph.get_sense(path) == handlegraph::PathSense::GENERIC) {
        // Generic paths only have a locus, so return whatever that is
        result.sample = graph.get_locus_name(path);
    } else {
        result.sample = graph.get_sample_name(path);
    }
    result.haplotype = graph.get_haplotype(path);

    return result;
}
}//end pangwas namespace

