#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <Eigen/Dense>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

#include "snarl_data_t.hpp"
#include "snarl_analyser.hpp"
#include "arg_parser.hpp"
#include "matrix.hpp"
#include "gaf_creator.hpp"
#include "post_processing.hpp"
namespace stoat_graph {

void print_help_graph();
int stoat_graph(int argc, char* argv[]);

}