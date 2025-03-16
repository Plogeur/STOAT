add_executable(tests binary_analysis_unit.cpp)
add_test(NAME MyTests COMMAND tests)

set(TEST_SOURCES
    binary_analysis_unit.cpp
    gaf_creator_unit.cpp
    list_snarl_paths_unit.cpp
    matrix_unit.cpp
    quantitative_analysis_unit.cpp
)

# Ensure test files can find the headers
include_directories(${CMAKE_SOURCE_DIR}/src)
