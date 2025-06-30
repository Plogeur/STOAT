#ifndef WRITER_INCLUDED
#define WRITER_INCLUDED

#include <iostream>
#include <handlegraph/path_position_handle_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>
#include "utils.hpp"

using namespace std;
namespace stoat_vcf{

    void write_binary_header(std::ofstream& outstream);
    void write_binary_covar_header(std::ofstream& outstream);
    void write_quantitative_header(std::ofstream& outstream);
    void write_eqtl_header(std::ofstream& outstream);


} //end namespace

#endif
