cmake_minimum_required(VERSION 3.10)
project(stoat VERSION 0.1.0 LANGUAGES CXX)

# C++ Standard Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "12.7") # IOS fix

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Find Packages
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# Copied from https://github.com/vgteam/vg/wiki/Programming-with-the-vg-API
# SDSL library ("sdsl" target) already exists from libbdsg
# GBWT doesn't have CMake yet so we need to explain how to build it.

# GBWT dependency: OpenMP
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # OpenMP may be installed by Homebrew in a place CMake can't find it by default.
    set (HOMEBREW_PREFIX "$ENV{HOMEBREW_PREFIX}")
    if ("${HOMEBREW_PREFIX}" STREQUAL "")
        # Ask Homebrew if it exists and where it is.
        execute_process(COMMAND brew --prefix OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE HOMEBREW_PREFIX)
    endif()
    if ("${HOMEBREW_PREFIX}" STREQUAL "")
        # Default to somewhere.
        set (HOMEBREW_PREFIX "/opt/homebrew")
    endif()
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}/opt/libomp")
endif()
find_package(OpenMP REQUIRED)

add_library(${PROJECT_NAME}_objs OBJECT
    src/utils.cpp
    src/arg_parser.cpp
    src/post_processing.cpp
    src/gaf_creator.cpp
    src/snarl_data_t.cpp
    src/matrix.cpp
    src/binary_table.cpp
    src/quantitative_table.cpp
    src/stats_test.cpp
    src/writer.cpp
    src/snarl_analyzer.cpp
    src/subcommand/vcf.cpp
    src/subcommand/graph.cpp
    #src/graph_path_association_finder.cpp
    #src/partitioner.cpp
)

# Include Directories
target_include_directories(${PROJECT_NAME}_objs PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
)

# Link Libraries
target_link_libraries(${PROJECT_NAME}_objs
    Boost::boost
    bdsg
    handlegraph
    Eigen3::Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# Add Executable
add_executable(${PROJECT_NAME} src/main.cpp)

# Link the code
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_objs)

target_compile_options(${PROJECT_NAME} PRIVATE -w) # fix warning

# Enable testing
enable_testing()

# Find Catch2 v3
find_package(Catch2 3 REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${Catch2_CMAKE_DIR})

include(CTest)
include(Catch)

# Add test executable
add_executable(unit_tests
    tests/unittest/binary_test_unit.cpp
    tests/unittest/gaf_creator_unit.cpp
    tests/unittest/snarl_data_t_unit.cpp
    tests/unittest/matrix_unit.cpp
    tests/unittest/quantitative_test_unit.cpp

    src/utils.cpp
    src/arg_parser.cpp
    src/gaf_creator.cpp
    src/snarl_data_t.cpp
    src/matrix.cpp
    src/binary_table.cpp
    src/quantitative_table.cpp
    src/stats_test.cpp
    src/writer.cpp
    src/snarl_analyzer.cpp
)

target_compile_options(unit_tests PRIVATE -w) # fix warning

# add src target lib + Catch2 lib 
target_link_libraries(unit_tests PRIVATE 
    Catch2::Catch2WithMain     
    Boost::boost
    bdsg
    handlegraph
    Eigen3::Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES})
