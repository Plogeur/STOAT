#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <Eigen/Dense>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

namespace stoat_graph {

void print_help_graph();
void stoat_graph(int argc, char* argv[]);

}