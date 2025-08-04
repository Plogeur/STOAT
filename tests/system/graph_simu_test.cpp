#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir) {
    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
    fs::create_directory(output_dir);
}

bool files_equal(const fs::path& f1, const fs::path& f2) {
    std::ifstream file1(f1, std::ios::binary);
    std::ifstream file2(f2, std::ios::binary);

    std::ostringstream buffer1, buffer2;
    buffer1 << file1.rdbuf();
    buffer2 << file2.rdbuf();

    return buffer1.str() == buffer2.str();
}

bool compare_output_dirs(const std::string& output_dir, const std::string& expected_dir) {
    for (const auto& file : fs::directory_iterator(expected_dir)) {
        auto expected_file = file.path();
        auto output_file = fs::path(output_dir) / expected_file.filename();

        if (!fs::exists(output_file)) {
            std::cerr << "Missing output file: " << output_file << "\n";
            return false;
        }
        if (!files_equal(expected_file, output_file)) {
            std::cerr << "Mismatch in file: " << expected_file.filename() << "\n";
            return false;
        }
    }
    return true;
}

bool run_test(
    const std::string& binary,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string phenotype_command,
    bool use_covariate = false) {

    clean_output_dir(output_dir);

    std::string cmd = binary + " graph"
        + " -g " + data_path + "/pg.full.pg"
        + " -d " + data_path + "/pg.full.dist"
        + " -S " + data_path + "/samples.g0.tsv"
        + phenotype_command + " -r ref";

    if (use_covariate) {
        cmd += " --covariate " + data_path + "/covariate.tsv"
             + " --covar-name CP1,SEX,CP3";
    }

    cmd += " --output " + output_dir;

    int result = std::system(cmd.c_str());
    if (result != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }

    return compare_output_dirs(output_dir, expected_dir);
}

TEST_CASE("Binary association tests graph", "[binary]") {
    const std::string binary = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/expected_output/graph/binary";
    const std::string expected_dir_covar = "../tests/expected_output/graph/binary_covar";
    const std::string data_path = "../data/binary";
    const std::string phenotype_command = " -T chi2 ";

    SECTION("Without covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, phenotype_command, false));
    }

    // SECTION("With covariate") {
    //     REQUIRE(run_test(binary, output_dir, expected_dir_covar, data_path, phenotype_command, true));
    // }
}

// TEST_CASE("Quantitative trait tests graph", "[quantitative]") {
//     const std::string binary = "./stoat";
//     const std::string output_dir = "../output_quantitative";
//     const std::string expected_dir = "../tests/expected_output/graph/quantitative";
//     const std::string expected_dir_covar = "../tests/expected_output/graph/quantitative_covar";
//     const std::string data_path = "../data/quantitative";
//     const std::string phenotype_command = " -q ";

//     //TODO: I added this so it would compile, idk what it should be
//     std::string sample_of_interest;

//     SECTION("Without covariate") {
//         REQUIRE(run_test(binary, output_dir, expected_dir, data_path, sample_of_interest, false));
//     }

//     SECTION("With covariate") {
//         REQUIRE(run_test(binary, output_dir, expected_dir_covar, data_path, sample_of_interest, true));
//     }
// }
