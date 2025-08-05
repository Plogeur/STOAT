#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

// === UTILITY FUNCTIONS ===

void clean_output_dir(const std::string& output_dir) {
    if (fs::exists(output_dir))
        fs::remove_all(output_dir);
    fs::create_directory(output_dir);
}

bool compare_output_files(const std::string& file1, const std::string& file2) {
    std::unordered_map<std::string, std::string> map1, map2;

    auto load_file = [](const std::string& path, std::unordered_map<std::string, std::string>& map) {
        std::ifstream infile(path);
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string token;
            std::vector<std::string> columns;

            while (std::getline(ss, token, '\t')) {
                columns.push_back(token);
            }

            if (columns.size() < 4) {
                throw std::runtime_error("Invalid line with fewer than 4 columns: " + line);
            }

            std::string snarl_id = columns[3]; // 4th column = SNARL
            map[snarl_id] = line;
        }
    };

    load_file(file1, map1);
    load_file(file2, map2);

    // Check missing in file2
    for (const auto& [snarl, line1] : map1) {
        auto it = map2.find(snarl);
        if (it == map2.end()) {
            std::cerr << "Missing SNARL in file2: " << snarl << std::endl;
            return false;
        } else if (line1 != it->second) {
            std::cerr << "Mismatch for SNARL " << snarl << ":\n"
                      << "File1: " << line1 << "\n"
                      << "File2: " << it->second << "\n";
            return false;
        }
    }

    // Check missing in file1
    for (const auto& [snarl, _] : map2) {
        if (map1.find(snarl) == map1.end()) {
            std::cerr << "Missing SNARL in file1: " << snarl << std::endl;
            return false;
        }
    }

    return true;
}

bool run_test_snarl(
    const std::string& binary,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string& phenotype) {

    clean_output_dir(output_dir);

    const std::string cmd = binary + " vcf"
        + " -p " + data_path + "/pg.full.pg"
        + " -d " + data_path + "/pg.full.dist"
        + " -r " + data_path + "/pg.chromosome"
        + " --output " + output_dir;

    int result = std::system(cmd.c_str());
    if (result != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return false;
    }
    
    std::string snarl_output = output_dir + '/' + phenotype + "_table.tsv";
    std::string snarl_expected = expected_dir + '/' + phenotype + "_table.tsv";

    return compare_output_files(snarl_output, snarl_expected);
}

bool run_test_gwas(
    const std::string& stoat_command,
    const std::string& output_dir,
    const std::string& expected_dir,
    const std::string& data_path,
    const std::string& phenotype,
    bool use_covariate = false) {

    std::string cmd = stoat_command + " vcf"
        + " -s " + output_dir + "/snarl_analyse.tsv"
        + " -v " + data_path + "/merged_output.vcf.gz";

    std::string type;

    if (phenotype == "eqtl") {
        cmd += " -e " + data_path + "/quantitative/qtl.tsv" 
        + " --gene-position " + data_path + "/quantitative/gene_position.tsv";

    } else if (phenotype == "binary") {
        cmd += " -b " + data_path + "/phenotype.tsv";

    } else if (phenotype == "quantitative") {
        cmd += " -q " + data_path + "/phenotype.tsv";
    }

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

    std::string gwas_output = output_dir + '/' + phenotype + "_table.tsv";
    std::string gwas_expected = expected_dir + '/' + phenotype + "_table.tsv";

    return compare_output_files(gwas_output, gwas_expected);
}

TEST_CASE("Binary association tests vcf", "[binary]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_binary";
    const std::string expected_dir = "../tests/expected_output/vcf/binary";
    const std::string expected_dir_covar = "../tests/expected_output/vcf/binary_covar";
    const std::string data_path = "../data/binary";
    const std::string phenotype = "binary";

    SECTION("Snarl decomposition") {
        REQUIRE(run_test_snarl(stoat_command, output_dir, expected_dir, data_path, phenotype));
    }

    SECTION("Without covariate") {
        REQUIRE(run_test_gwas(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test_gwas(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

TEST_CASE("Quantitative trait tests vcf", "[quantitative]") {
    const std::string stoat_command = "../bin/stoat";
    const std::string output_dir = "../output_quantitative";
    const std::string expected_dir = "../tests/expected_output/vcf/quantitative";
    const std::string expected_dir_covar = "../tests/expected_output/vcf/quantitative_covar";
    const std::string data_path = "../data/quantitative";
    const std::string phenotype = "quantitative";

    SECTION("Snarl decomposition") {
        REQUIRE(run_test_snarl(stoat_command, output_dir, expected_dir, data_path, phenotype));
    }

    SECTION("Without covariate") {
        REQUIRE(run_test_gwas(stoat_command, output_dir, expected_dir, data_path, phenotype, false));
    }

    SECTION("With covariate") {
        REQUIRE(run_test_gwas(stoat_command, output_dir, expected_dir_covar, data_path, phenotype, true));
    }
}

// TEST_CASE("eQTL tests vcf", "[eqtl]") {
//     const std::string stoat_command = "./stoat";
//     const std::string output_dir = "../output_eqtl";
//     const std::string expected_dir = "../tests/expected_output/vcf/eqtl";
//     const std::string expected_dir_covar = "../tests/expected_output/vcf/eqtl_covar";
//     const std::string phenotype_command = " -e ";

//     const std::string data_path = "../data/eqtl";

//     SECTION("Snarl decomposition") {
//         REQUIRE(run_test_snarl(stoat_command, output_dir, expected_dir, data_path, phenotype));
//     }

//     SECTION("Without covariate") {
//         REQUIRE(run_test_gwas(stoat_command, output_dir, expected_dir, data_path, phenotype_command, false, true));
//     }

//     SECTION("With covariate") {
//         REQUIRE(run_test_gwas(stoat_command, output_dir, expected_dir_covar, data_path, phenotype_command, true, true));
//     }
// }
