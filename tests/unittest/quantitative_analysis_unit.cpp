#define CATCH_CONFIG_MAIN unit_test::quantitative_analysis_unit
#include "catch.hpp"
#include "quantitative_analysis.hpp"

TEST_CASE("Linear Regression Test", "[linear_regression]") {
    // Create a sample dataframe and quantitative phenotype
    std::unordered_map<std::string, std::vector<int>> df = {
        {"Sample1", {1, 2, 3}},
        {"Sample2", {4, 5, 6}},
        {"Sample3", {7, 8, 9}}
    };

    std::unordered_map<std::string, double> quantitative_phenotype = {
        {"Sample1", 10.0},
        {"Sample2", 20.0},
        {"Sample3", 30.0}
    };

    // Run the linear regression function
    auto [se, beta, p_value] = linear_regression(df, quantitative_phenotype);

    // Check the returned results for expected values
    REQUIRE(se != 0.0);
    REQUIRE(beta != 0.0);
    REQUIRE(p_value != "NA");
    REQUIRE(p_value != "0.0000");  // Adjust this to an expected p-value if you know the data
}

// TEST_CASE("Quantitative Table Creation Test", "[create_quantitative_table]") {
//     // Create dummy data
//     std::vector<std::string> list_samples = {"Sample1", "Sample2", "Sample3"};
//     std::vector<std::string> column_headers = {"Path1", "Path2"};
    
//     // Mock Matrix class or use a proper one if available
//     Matrix matrix; // Define a mock matrix with the required interface
    
//     // Create the quantitative table
//     std::unordered_map<std::string, std::vector<int>> df = create_quantitative_table(list_samples, column_headers, matrix);
    
//     // Verify that the output dataframe is correctly populated
//     REQUIRE(df.size() == list_samples.size());
    
//     // Verify that each sample in the dataframe has the expected number of columns (genotype counts)
//     for (const auto& sample : df) {
//         REQUIRE(sample.second.size() == column_headers.size());
//     }
// }
