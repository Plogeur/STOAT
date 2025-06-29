#include "association_finder.hpp"
#include "utils.hpp"

//#define DEBUG_ASSOCIATION_FINDER

namespace pangwas {

AssociationFinder::AssociationFinder(const handlegraph::PathPositionHandleGraph& graph, 
                                     const bdsg::SnarlDistanceIndex& distance_index,
                                     std::string test_method,
                                     const std::set<std::string>& samples_of_interest, 
                                     std::string reference_name, 
                                     std::string output_format, 
                                     std::ostream& out_associated, 
                                     std::ostream& out_unassociated, 
                                     size_t allele_size_limit, double p_value) :
    graph(graph), 
    distance_index(distance_index), 
    samples_of_interest(samples_of_interest), 
    reference_name(reference_name), 
    output_format(output_format),
    out_associated(out_associated),
    out_unassociated(out_unassociated),
    allele_size_limit(allele_size_limit)
    {}

std::pair<bool, std::unordered_set<std::string>> AssociationFinder::is_snarl_associated(const handlegraph::net_handle_t& snarl) const {

    std::vector<std::set<std::string>> sample_partitions = partition_samples_in_snarl(snarl);


    bool is_associated = false;
    std::unordered_set<std::string> samples_to_print;
    for (std::set<std::string>& partition : sample_partitions) {
        if (tester->is_associated(partition)) {
            is_associated = true;
            samples_to_print.insert(std::make_move_iterator(partition.begin()), std::make_move_iterator(partition.end()));
        } else {
            samples_to_print.insert(std::make_move_iterator(partition.begin()), std::make_move_iterator(++partition.begin()));
        }
    }
    return std::make_pair(is_associated, std::move(samples_to_print));
}


void AssociationFinder::write_associated_snarls() const {

    // If the file output has a header, write it
    write_header();


    std::vector<handlegraph::net_handle_t> chains;
    chains.reserve(graph.get_node_count()/100);
    handlegraph::net_handle_t root = distance_index.get_root();
    distance_index.for_each_child(root, [&] (handlegraph::net_handle_t chain) {
        chains.emplace_back(chain);
        return true;
    });

    while (!chains.empty()) {
        handlegraph::net_handle_t chain = chains.back();
        chains.pop_back();

        distance_index.for_each_child(chain, [&] (handlegraph::net_handle_t snarl) {

            //TODO: For now it's fine to check is_eligible here because it's only checking size and we don't want to look at small chains anyway
            if (distance_index.is_snarl(snarl) && snarl_is_eligible(snarl) ) {
                std::pair<bool, std::unordered_set<std::string>> associated = is_snarl_associated(snarl);
                if (associated.first) {
                    // If the snarl is associated, output it

                    write_snarl(snarl, associated.second);

                } else {

                    // If it is not associated, add its child chains to the stack
                    distance_index.for_each_child(snarl, [&] (handlegraph::net_handle_t child) {

                        chains.emplace_back(child);
                        return true;
                    });
                }
            }
            return true;
        });
    }
}

bool AssociationFinder::snarl_is_eligible(const handlegraph::net_handle_t& snarl) const {
    return distance_index.maximum_length(snarl) >= allele_size_limit;
}

void AssociationFinder::write_header() const {
    if (output_format == "tsv") {
        out_associated << "path_name\tstart_offset\tend_offset\tvariant_size" << std::endl;
    }
    //TODO: When adding stuff, check if out_associated/unassociated are the same and write to both if necessary
}

void AssociationFinder::write_snarl(const handlegraph::net_handle_t& snarl, const std::unordered_set<std::string>& samples) const {
    if (output_format == "tsv") {
        write_tsv_of_snarl(snarl);
    } else if (output_format == "fasta") {
        //Write one fasta record for each path set
        write_fasta_of_snarl(snarl, samples);
    } else {
       std::cerr << "error[pangwas]: unknown output format " << output_format << std::endl;
    }
}

void AssociationFinder::write_tsv_of_snarl(const handlegraph::net_handle_t& snarl) const {
    // Write all coordinates for a single path
    handlegraph::path_handle_t ref_path;
    bool first = true;
    for (const path_range_t& range : get_coordinates_of_snarl(snarl, true, reference_name, false)) {
        if (first) {
            first = false;
            ref_path = graph.get_path_handle_of_step(range.start);
        } else if (graph.get_path_handle_of_step(range.start) != ref_path) {
            continue;
        }

        out_associated << graph.get_path_name(graph.get_path_handle_of_step(range.start)) << "\t"
            << (graph.get_position_of_step(range.start) + 
                distance_index.minimum_length(distance_index.get_net(graph.get_handle_of_step(range.start), &graph))) << "\t"
            << graph.get_position_of_step(range.end) << "\t"
            << distance_index.maximum_length(snarl) 
            << std::endl;
    }
}


void AssociationFinder::write_fasta_of_snarl(const handlegraph::net_handle_t& snarl, const std::unordered_set<std::string>& samples) const {

    // A handle_t of the start bound facing in
    handlegraph::handle_t start_handle = distance_index.get_handle(distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true)), &graph);

    // Get a unique name for the snarl, as the start and end ids
    // I think even in the case of a looping chain where there is another snarl on the other end, the order of node ids will be flipped
    std::string snarl_name = "snarl:" + 
                             std::to_string((int)graph.get_id(start_handle)) + 
                             "-" + 
                             std::to_string((int)distance_index.node_id(distance_index.get_bound(snarl, true, false)));

    // Get a reference range for the snarl.
    // If the reference goes through the snarl multiple times, get the largest interval

    std::vector<path_range_t> ref_ranges = get_coordinates_of_snarl(snarl, true, reference_name, false);
    std::string ref_coordinates  = "NOREF:?:?";
    int start_offset = std::numeric_limits<int>::max();
    int end_offset = 0;
    //Only get the coordinates for one path
    bool first = true;
    handlegraph::path_handle_t ref_path;
    for (const path_range_t& ref_range : ref_ranges){
        if (first) {
            first = false;
            ref_path = graph.get_path_handle_of_step(ref_range.start);
        } else if (graph.get_path_handle_of_step(ref_range.start) != ref_path) {
            continue;
        }
        ref_coordinates = graph.get_path_name(graph.get_path_handle_of_step(ref_range.start));
        start_offset = std::min(start_offset, 
                                (int)(graph.get_position_of_step(ref_range.start) + 
                                      distance_index.minimum_length(distance_index.get_net(graph.get_handle_of_step(ref_range.start), &graph))));
        end_offset = std::max(end_offset, (int)graph.get_position_of_step(ref_range.end));
    }

    if (ref_ranges.size() != 0) {
        ref_coordinates += ":" + std::to_string(start_offset) + "-" + std::to_string(end_offset);
    }

    // Now go through each path that goes through the snarl and print the sequence
    std::vector<path_range_t> path_ranges = get_coordinates_of_snarl(snarl, false, "", true);
    for (const path_range_t& path_range : path_ranges) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(path_range.start);
        if (samples.empty() || samples.count(get_sample_name_from_path(graph, path)) != 0) {
            //If we aren't checking samples, or if this is a sample we want

            // Print to the correct stream, depending on if it is associated or not
            std::ostream& out = samples_of_interest.count(get_sample_name_from_path(graph, path)) != 0
                              ? out_associated
                              : out_unassociated;

            // Print the header
            out << ">" << snarl_name << "|" 
                << ref_coordinates << "|" 
                << graph.get_path_name(path) << ":" 
                << (graph.get_position_of_step(path_range.start) + distance_index.minimum_length(distance_index.get_net(graph.get_handle_of_step(path_range.start), &graph))) << "-" 
                << graph.get_position_of_step(path_range.end) << std::endl; 

            // Now print the sequence in 80bp chunks.
            // Keep a buffer to print 80 bp at a time
            std::string sequence_buffer = ""; 
            handlegraph::step_handle_t next_step = graph.get_next_step(path_range.start);
            while (next_step != path_range.end) {
                std::string node_seq = graph.get_sequence(graph.get_handle_of_step(next_step));
                while (node_seq.size() != 0) {

                    // Fill in sequence_buffer up to 80 characters
                    size_t to_add = 80 - sequence_buffer.size();
                    sequence_buffer += node_seq.substr(0, to_add);
                    node_seq.erase(0, to_add);

                    // If the buffer is full, write it and clear it
                    if (sequence_buffer.size() == 80) {
                        out << sequence_buffer << std::endl;
                        sequence_buffer.clear();
                    }
                }
                handlegraph::step_handle_t step = next_step;
                if (!graph.has_next_step(step)) {
                    break;
                }
                next_step = graph.get_next_step(step);
            }
            if (!sequence_buffer.empty()) {
                out << sequence_buffer << std::endl;
            }
        }
    }
}

std::vector<AssociationFinder::path_range_t> AssociationFinder::get_coordinates_of_snarl(const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths) const {
    std::vector<AssociationFinder::path_range_t> ranges;
    // If a sample name is given, then always look for that first
    if (!sample_name.empty()) {
        ranges = get_coordinates_of_snarl_helper(snarl, false, sample_name, false);
        if (!ranges.empty()) {
            return ranges;
        }
    }
    if (get_reference) {
        //If we didn't find the specific path and we are looking for a reference, look for a reference-sense path next

        //Try with reference-sense path
        ranges = get_coordinates_of_snarl_helper(snarl, true, "", false);
        if (!ranges.empty()) {
            return ranges;
        }
    }
    if (get_all_paths) {
        //If we just want all paths, return that 

        ranges = get_coordinates_of_snarl_helper(snarl, false, "", true);
        return ranges;

    } else {
        //Try with any path

        handlegraph::net_handle_t start_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true));
        handlegraph::net_handle_t end_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, true));

        //Get all paths on the start node
        std::set<handlegraph::path_handle_t> start_paths;
        graph.for_each_step_on_handle(distance_index.get_handle(start_net, &graph), [&](const handlegraph::step_handle_t& step) {
            start_paths.insert(graph.get_path_handle_of_step(step));
            return true;
        });

        std::string new_sample_name = "";
        //Now go through the end node and find a path that also goes through the start node
        graph.for_each_step_on_handle(distance_index.get_handle(end_net, &graph), [&](const handlegraph::step_handle_t& step) {
            if (start_paths.count(graph.get_path_handle_of_step(step)) != 0) {
                new_sample_name = graph.get_path_name(graph.get_path_handle_of_step(step));
                return false;
            }
            return true;
        });

        if (new_sample_name.empty()) {
            return ranges;
        } else {
            ranges = get_coordinates_of_snarl_helper(snarl, false, new_sample_name, false);
            return ranges;
        }

    }
}
std::vector<AssociationFinder::path_range_t> AssociationFinder::get_coordinates_of_snarl_helper(const handlegraph::net_handle_t& snarl, bool get_reference, std::string sample_name, bool get_all_paths) const {
    #ifdef DEBUG_ASSOCIATION_FINDER
   std::cerr << "Get coordinates of " << distance_index.net_handle_as_string(snarl) << std::endl;
    if (get_reference) {
        assert(sample_name.empty());
        assert(!get_all_paths);
    }
    if (!sample_name.empty()) {
        assert(!get_reference);
        assert(!get_all_paths);
    }
    if (get_all_paths) {
        assert(!get_reference);
        assert(sample_name.empty());
    }
    #endif

    // Bound nodes going into of the snarl
    handlegraph::net_handle_t start_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, false, true));
    handlegraph::net_handle_t end_net = distance_index.get_node_from_sentinel(distance_index.get_bound(snarl, true, true));

    // Map path to the steps on the path that traverse the snarl bounds
    std::map<handlegraph::path_handle_t, std::vector<handlegraph::step_handle_t>> path_to_steps;

    // Keep track if we found a traversal of the snarl (may be start-start or end-end)
    bool found_pair = false;

    // Get the step_handles, filtering for the paths we are interested in
    // Steps don't care about the orientation of the handle, they will always (I think) be going forwards in the path
    graph.for_each_step_on_handle(distance_index.get_handle(start_net, &graph), [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);
        if ((get_reference && (graph.get_sense(path) == handlegraph::PathSense::REFERENCE)) ||
            (!sample_name.empty() && graph.get_path_name(path).find(sample_name) != std::string::npos) ||
            (get_all_paths)) {
            // If we are looking for a reference path and this is a reference path
            // or if this is the path we want or if we want all paths
            if (path_to_steps.count(path) == 0) {
                path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
            } else {
                found_pair = true;
            }
            path_to_steps[path].emplace_back(step);
        }
        return true;
    });
    #ifdef DEBUG_ASSOCIATION_FINDER
   std::cerr << "After start node, found" << std::endl;
    for (const auto& x : path_to_steps) {
       std::cerr << graph.get_path_name(x.first) << ": " << x.second.size() << std::endl;
    }
    #endif
    graph.for_each_step_on_handle(distance_index.get_handle(end_net, &graph), [&] (const handlegraph::step_handle_t& step) {
        handlegraph::path_handle_t path = graph.get_path_handle_of_step(step);
        if ((get_reference && (graph.get_sense(path) == handlegraph::PathSense::REFERENCE)) ||
            (!sample_name.empty() && graph.get_path_name(path).find(sample_name) != std::string::npos) ||
            (get_all_paths)) {
            // If we are looking for a reference path and this is a reference path
            // or if this is the path we are looking for, or we want all paths 
            if (path_to_steps.count(path) == 0) {
                path_to_steps[path] = std::vector<handlegraph::step_handle_t>();
            } else {
                found_pair = true;
            }
            path_to_steps[path].emplace_back(step);
        }
        return true;
    });

    #ifdef DEBUG_ASSOCIATION_FINDER
   std::cerr << "After end node, found" << std::endl;
    for (const auto& x : path_to_steps) {
       std::cerr << graph.get_path_name(x.first) << ": " << x.second.size() << std::endl;
    }
    #endif

    std::vector<path_range_t> ranges;

    if (found_pair) {
        //If we found a path going through the snarl, return the pairs

        for (auto& path_steps : path_to_steps) {
            const handlegraph::path_handle_t& path = path_steps.first;
            std::vector<handlegraph::step_handle_t>& steps = path_steps.second;

            if (steps.size() < 2) {
                continue;
            }
            std::sort(steps.begin(), steps.end(), [&] (const handlegraph::step_handle_t& a, const handlegraph::step_handle_t& b) {
                return graph.get_position_of_step(a) < graph.get_position_of_step(b);
            });

            #ifdef DEBUG_ASSOCIATION_FINDER
            for (size_t step_i = 0 ; step_i < steps.size() ; step_i++) {
                if (step_i % 2 == 0) {
                    // If this is an even number, then the path should go into the snarl
                    assert(graph.get_handle_of_step(steps[step_i]) == distance_index.get_handle(start_net, &graph) ||
                        graph.get_handle_of_step(steps[step_i]) == distance_index.get_handle(end_net, &graph));
                } else {
                    //If this is an odd number, it should go out of the snarl
                    assert(graph.get_handle_of_step(steps[step_i]) == graph.flip(distance_index.get_handle(start_net, &graph)) ||
                        graph.get_handle_of_step(steps[step_i]) == graph.flip(distance_index.get_handle(end_net, &graph)));
                }
            }
            #endif
            for (size_t i = 0 ; i < steps.size() ; i += 2) {
                ranges.push_back({steps[i], steps[i+1]});
            }
        }
        return ranges;

    } else {
        // If we didn't find anything
        return ranges;
    }

}




}//end pangwas namespace

