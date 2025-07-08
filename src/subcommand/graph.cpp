// This file is part of STOAT 0.0.1, copyright (C) 2024-2025 
// Authors : Matis Alias-Bagarre, Xian-hui Chang & Jean Monlong.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <cstdlib>
#include <getopt.h>
#include <omp.h>

#include "../snarl_data_t.hpp"
#include "../snarl_analyser.hpp"
#include "../arg_parser.hpp"
#include "../matrix.hpp"
#include "../gaf_creator.hpp"
#include "../post_processing.hpp"

namespace stoat_command {

const static void print_help_graph() {
    std::cerr   << "usage: stoat graph [options]\n\n"
                << "  -g, --graph FILE                   use this graph (only hash graph works for now) (required)" << endl
                << "  -d, --distance-index FILE          use this distance index (required)" << endl
                << "  -s, --sample-of-interest NAME      the name of the sample with the trait of interest (may repeat)" << endl
                << "  -o, --output-format NAME           the format of the output (tsv / fasta) [tsv]" << endl
                << "  -a, --associated-filename FILE     write the records for the associated samples to FILE" << endl
                << "  -u, --unassociated-filename FILE   write the records for the unassociated samples to FILE" << endl
                << "  -t, --threads N                    number of threads to use" << endl
                << "  -T, --test NAME                    which test will be used to determine association (exact / fishers / chi2) [exact]" << endl
                << "  -p, --p-value-threshold FLOAT      what is the threshold p-value to be considered significant? [0.05]" << endl
                << "  -m, --method NAME                  what method is used to find associations? (paths) [paths]" << endl
                << "  -l, --allele-size-limit INT        don't report variants smaller than this [0]" << endl
                << "  -r, --reference-sample NAME        if there is no reference in the graph, use this sample as the reference" << endl
                << "  -h, --help                         print this help message" << std::endl;
}

const static int main_graph(int argc, char* argv[]) {

    std::string graph_name;
    std::string distance_name;
    size_t allele_size_limit = 0;
    double p_value = 0.05;
    std::string method_name = "paths";
    std::string test_method = "exact";
    std::string reference_sample;
    std::set<std::string> samples_of_interest;
    std::string output_format= "tsv";
    std::string associated_filename;
    std::string unassociated_filename;

//     int c = 0;
//     optind = 1;
//     while (true) {
//         static struct option long_options[] =
//             {
//                 {"graph", required_argument, 0, 'g'},
//                 {"distance-index", required_argument, 0, 'd'},
//                 {"allele-size-limit", required_argument, 0, 'l'},
//                 {"threads", required_argument, 0, 't'},
//                 {"test", required_argument, 0, 'T'},
//                 {"p-value", required_argument, 0, 'p'},
//                 {"method", required_argument, 0, 'm'},
//                 {"reference-sample", required_argument, 0, 'r'},
//                 {"sample-of-interest", required_argument, 0, 's'},
//                 {"output-format", required_argument, 0, 'o'},
//                 {"associated-filename", required_argument, 0, 'a'},
//                 {"unassociated-filename", required_argument, 0, 'u'},
//                 {"help", no_argument, 0, 'h'},
//                 {0, 0, 0, 0}
//             };

//         int option_index = 0;
//         c = getopt_long(argc, argv, "g:d:l:t:T:p:m:r:s:o:a:u:h",
//                         long_options, &option_index); 
//         if (c == -1) {
//             break;
//         }
//         switch (c) {
//             case 'g':
//                 graph_name = optarg;
//                 break;
//             case 'd':
//                 distance_name = optarg;
//                 break;
//             case 'l':
//                 allele_size_limit = std::stoi(optarg);
//                 break;
//             case 't':
//                 omp_set_num_threads(std::stoi(optarg));
//                 break;
//             case 'T':
//                 test_method = optarg;
//                 break;
//             case 'p':
//                 p_value = std::stof(optarg);
//                 break;
//             case 'm':
//                 method_name = optarg;
//                 break;
//             case 'r':
//                 reference_sample = optarg;
//                 break;
//             case 's':
//                 samples_of_interest.emplace(optarg);
//                 break;
//             case 'o':
//                 output_format = optarg;
//                 break;
//             case 'a':
//                 associated_filename = optarg;
//                 break;
//             case 'u':
//                 unassociated_filename = optarg;
//                 break;
//             case 'h': print_help_graph(); exit(EXIT_SUCCESS); break;
//             default:
//                 std::cerr << "Unknown argument. Use -h or --help for usage.\n";
//                 return EXIT_FAILURE;
//         }
//     }

//     // Check that the inputs are ok
//     if (graph_name.empty()) {
//         std::cerr << "error [pangwas]: pangwas requires a graph file" << std::endl;
//         EXIT_FAILURE; 
//     }
//     if (distance_name.empty()) {
//         std::cerr << "error [pangwas]: pangwas requires a distance index file" << std::endl;
//         EXIT_FAILURE; 
//     }
//     if (samples_of_interest.empty()) {
//         std::cerr << "error [pangwas]: pangwas requires samples of interest" << std::endl;
//         EXIT_FAILURE; 
//     }

//     // Tell the IO library about libvg types.
//     if (!pangwas::io::register_libvg_io()) {
//        std::cerr << "error[vg]: Could not register libvg types with libvgio" << std::endl;
//         EXIT_FAILURE;
//     }

//     // Load the graph and make it a PathPositionHandleGraph
//     unique_ptr<handlegraph::PathHandleGraph> path_graph = vg::io::VPKG::load_one<handlegraph::PathHandleGraph>(graph_name);
//     bdsg::PathPositionOverlayHelper overlay_helper;
//     bdsg::PathPositionHandleGraph* graph = overlay_helper.apply(path_graph.get());

//     // Load the distance index
//     bdsg::SnarlDistanceIndex distance_index;
//     distance_index.deserialize(distance_name);

//     // Get the out streams
//     std::ofstream out_associated;
//     if (!associated_filename.empty()) {
//         out_associated.open(associated_filename);
//     }
//     std::ofstream out_unassociated;
//     if (!unassociated_filename.empty()) {
//         out_unassociated.open(unassociated_filename);
//     }

//     if (method_name == "paths") {
//         pangwas::PathAssociationFinder af (*graph, 
//                                         distance_index, 
//                                         test_method,
//                                         samples_of_interest, 
//                                         reference_sample, 
//                                         output_format,
//                                         associated_filename.empty() ? std::cout : out_associated,
//                                         unassociated_filename.empty() ? std::cout : out_unassociated,
//                                         allele_size_limit,
//                                         p_value);
//         af.write_associated_snarls();
//     } else {
//         std::cerr << "error [pangwas]: unknown method " << method_name << std::endl;
//         EXIT_FAILURE; 
//     }

//     //Close streams
//     if (!associated_filename.empty()) {
//         out_associated.close();
//     }
//     if (!unassociated_filename.empty()) {
//         out_unassociated.close();
//     }

    return EXIT_SUCCESS;
}

}
