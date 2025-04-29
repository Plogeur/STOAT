#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <unordered_map>
#include <string>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <tuple>
#include <vector>
#include "../../src/arg_parser.hpp"  // adjust path as needed

// Helper to create a temporary test file
std::string create_test_pheno_file(const std::string& content, const std::string& filename = "test_pheno.txt") {
    std::ofstream out(filename);
    out << content;
    out.close();
    return filename;
}

// TEST_CASE("Binary phenotype parsing", "[parse_binary_pheno]") {
//     SECTION("Valid phenotype file with mixed cases and controls") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1 1\n"
//             "F2 I2 2\n"
//             "F3 I3 1\n"
//             "F4 I4 2\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         auto result = parse_binary_pheno(file_path);

//         REQUIRE(result.size() == 4);
//         REQUIRE(result["I1"] == false);
//         REQUIRE(result["I2"] == true);
//         REQUIRE(result["I3"] == false);
//         REQUIRE(result["I4"] == true);
//     }

//     SECTION("Invalid header") {
//         std::string file_content =
//             "FID ID PHENO\n"
//             "F1 I1 1\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_binary_pheno(file_path), std::invalid_argument);
//     }

//     SECTION("Non-binary phenotype value") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1 3\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_binary_pheno(file_path), std::runtime_error);
//     }

//     SECTION("Malformed line") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_binary_pheno(file_path), std::runtime_error);
//     }

//     SECTION("Non-integer phenotype") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1 X\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_binary_pheno(file_path), std::runtime_error);
//     }
// }


// TEST_CASE("Quantitative phenotype parsing", "[parse_quantitative_pheno]") {
//     SECTION("Valid quantitative phenotype file") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1 1.5\n"
//             "F2 I2 2.0\n"
//             "F3 I3 -3.25\n"
//             "F4 I4 0.0\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         auto result = parse_quantitative_pheno(file_path);

//         REQUIRE(result.size() == 4);
//         REQUIRE(result["I1"] == 1.5);
//         REQUIRE(result["I2"] == 2.0);
//         REQUIRE(result["I3"] == -3.25);
//         REQUIRE(result["I4"] == 0.0);
//     }

//     SECTION("Invalid header in quantitative phenotype file") {
//         std::string file_content =
//             "FID ID PHENO\n"
//             "F1 I1 1.5\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_quantitative_pheno(file_path), std::invalid_argument);
//     }

//     SECTION("Non-numeric phenotype value") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1 abc\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_quantitative_pheno(file_path), std::runtime_error);
//     }

//     SECTION("Malformed line in quantitative phenotype file") {
//         std::string file_content =
//             "FID IID PHENO\n"
//             "F1 I1\n";
//         std::string file_path = create_test_pheno_file(file_content);

//         REQUIRE_THROWS_AS(parse_quantitative_pheno(file_path), std::runtime_error);
//     }
// }

// // Helper to write minimal VCF content to file
// void write_vcf_file(const std::string& path, const std::string& content) {
//     std::filesystem::create_directories("test_data");
//     std::ofstream file(path);
//     file << content;
//     file.close();
// }

// TEST_CASE("VCF Parsing", "[parse_vcf][parseHeader]") {
//     SECTION("Valid VCF with 3 samples") {
//         std::string vcf_content =
//             "##fileformat=VCFv4.2\n"
//             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
//             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsampleA\tsampleB\tsampleC\n"
//             "1\t12345\trsTest\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n";

//         std::string vcf_path = "test_data/test_valid.vcf";
//         write_vcf_file(vcf_path, vcf_content);

//         auto [samples, ptr, hdr, rec] = parseHeader(vcf_path);

//         REQUIRE(samples.size() == 3);
//         REQUIRE(samples[0] == "sampleA");
//         REQUIRE(samples[1] == "sampleB");
//         REQUIRE(samples[2] == "sampleC");

//         bcf_destroy(rec);
//         bcf_hdr_destroy(hdr);
//         bcf_close(ptr);
//     }

//     SECTION("Malformed VCF should throw error") {
//         std::string vcf_content =
//             "This is not a VCF\n"
//             "Just some junk lines\n";

//         std::string vcf_path = "test_data/test_invalid.vcf";
//         write_vcf_file(vcf_path, vcf_content);

//         REQUIRE_THROWS_AS(parseHeader(vcf_path), std::runtime_error);
//     }

//     SECTION("Empty VCF should throw error on header read") {
//         std::string vcf_path = "test_data/test_empty.vcf";
//         write_vcf_file(vcf_path, "");  // Empty file

//         REQUIRE_THROWS_AS(parseHeader(vcf_path), std::runtime_error);
//     }
// }

// TEST_CASE("Check sample matching", "[check_match_samples]") {
//     SECTION("All keys found, size match (bool)") {
//         std::unordered_map<std::string, bool> pheno = {
//             {"sample1", true},
//             {"sample2", false},
//             {"sample3", true}
//         };
//         std::vector<std::string> vcf_samples = {"sample1", "sample2", "sample3"};

//         REQUIRE_NOTHROW(check_match_samples<bool>(pheno, vcf_samples));
//     }

//     SECTION("All keys found, size mismatch (double)") {
//         std::unordered_map<std::string, double> pheno = {
//             {"sample1", 1.0},
//             {"sample2", 2.0},
//             {"sample3", 3.0},
//             {"extra_sample", 4.0}
//         };
//         std::vector<std::string> vcf_samples = {"sample1", "sample2", "sample3"};

//         // Should not throw, just a warning (which we can't easily assert here)
//         REQUIRE_NOTHROW(check_match_samples<double>(pheno, vcf_samples));
//     }

//     SECTION("Missing key throws runtime_error") {
//         std::unordered_map<std::string, bool> pheno = {
//             {"sample1", true},
//             {"sample2", false}
//         };
//         std::vector<std::string> vcf_samples = {"sample1", "sample2", "sample3"};

//         REQUIRE_THROWS_AS(check_match_samples<bool>(pheno, vcf_samples), std::runtime_error);
//     }
// }

// // Helper function to create a test file
// std::string create_test_covar_file(const std::string& content, const std::string& filename = "test_covariate.txt") {
//     std::ofstream out(filename);
//     out << content;
//     out.close();
//     return filename;
// }

// TEST_CASE("Covariate format checking", "[check_format_covariate]") {
//     SECTION("Consistent column count") {
//         std::string content =
//             "IID age sex\n"
//             "A 25 1\n"
//             "B 30 0\n"
//             "C 45 1\n";
//         std::string path = create_test_covar_file(content);

//         REQUIRE_NOTHROW(check_format_covariate(path));
//     }

//     SECTION("Inconsistent column count prints error") {
//         std::string content =
//             "IID age sex\n"
//             "A 25 1\n"
//             "B 30\n" // Missing one column
//             "C 45 1\n";
//         std::string path = create_test_covar_file(content);
//         REQUIRE_THROWS_AS(check_format_covariate(path), std::runtime_error);
//     }
// }

// TEST_CASE("Parse covariates from file", "[parse_covariates]") {
//     SECTION("Valid covariate file and columns") {
//         std::string content =
//             "IID age sex\n"
//             "A 25 1\n"
//             "B 30 0\n"
//             "C 45 1\n";
//         std::string path = create_test_covar_file(content);
//         std::vector<std::string> covars = {"age", "sex"};

//         auto result = parse_covariates(path, covars);

//         REQUIRE(result.size() == 3);
//         REQUIRE(result["A"][0] == 25);
//         REQUIRE(result["A"][1] == 1);
//         REQUIRE(result["B"][0] == 30);
//         REQUIRE(result["B"][1] == 0);
//         REQUIRE(result["C"][0] == 45);
//     }

//     SECTION("Missing IID column") {
//         std::string content =
//             "ID age sex\n"
//             "A 25 1\n";
//         std::string path = create_test_covar_file(content);
//         std::vector<std::string> column_covars = {"age", "sex"};
//         REQUIRE_THROWS_AS(parse_covariates(path, column_covars), std::runtime_error);
//     }

//     SECTION("Missing covariate column") {
//         std::string content =
//             "IID age sex\n"
//             "A 25 1\n";
//         std::string path = create_test_covar_file(content);
//         std::vector<std::string> column_covars = {"height"}; // not present
//         REQUIRE_THROWS_AS(parse_covariates(path, column_covars), std::runtime_error);
//     }

//     SECTION("Non-numeric value in covariate field") {
//         std::string content =
//             "IID age sex\n"
//             "A XX 1\n"; // XX is not numeric
//         std::string path = create_test_covar_file(content);
//         std::vector<std::string> column_covars = {"age", "sex"};
//         REQUIRE_THROWS_AS(parse_covariates(path, column_covars), std::runtime_error);
//     }
// }
