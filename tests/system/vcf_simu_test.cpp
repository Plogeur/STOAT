#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

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
    bool use_covariate = false) {

    clean_output_dir(output_dir);

    std::string cmd = binary + " vcf"
        + " -p " + data_path + "/pg.pg"
        + " -d " + data_path + "/pg.dist"
        + " -v " + data_path + "/merged_output.vcf.gz"
        + " -b " + data_path + "/phenotype.tsv";

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

TEST_CASE("Binary association tests", "[binary]") {
    const std::string binary = "./stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../vcf/expected_output/binary";
    const std::string data_path = "../data/binary";

    SECTION("Without covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, true));
    }
}

TEST_CASE("Quantitative trait tests", "[quantitative]") {
    const std::string binary = "./stoat";
    const std::string output_dir = "../output_quantitative";
    const std::string expected_dir = "../vcf/expected_output/quantitative";
    const std::string data_path = "../data/quantitative";

    SECTION("Without covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, true));
    }
}

TEST_CASE("eQTL tests", "[eqtl]") {
    const std::string binary = "./stoat";
    const std::string output_dir = "../output_eqtl";
    const std::string expected_dir = "../vcf/expected_output/eqtl";
    const std::string data_path = "../data/eqtl";

    SECTION("Without covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test(binary, output_dir, expected_dir, data_path, true));
    }
}
