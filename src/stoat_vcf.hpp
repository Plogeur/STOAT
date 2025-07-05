#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <Eigen/Dense>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

namespace stoat_vcf {

void print_help_vcf();
void stoat_vcf(int argc, char* argv[]);

}