cmake_minimum_required(VERSION 3.10)
project(stoat VERSION 0.1.1 LANGUAGES CXX)

# C++ Standard Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "12.7") # IOS fix

include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Find Packages
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# Copied from https://github.com/vgteam/vg/wiki/Programming-with-the-vg-API
# SDSL library ("sdsl" target) already exists from libbdsg

# Detect Homebrew libomp
execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)

# Force OpenMP Flags for Apple Clang
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")

find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP")
    add_compile_options(${OpenMP_CXX_FLAGS})
    link_directories(${LIBOMP_PREFIX}/lib)
    link_libraries(${OpenMP_omp_LIBRARY})
else()
    message(FATAL_ERROR "OpenMP not found or misconfigured.")
endif()

# Add Executable
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/stoat_vcf.cpp
    src/stoat_graph.cpp
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
    src/snarl_analyser.cpp
    ${ALGLIB_SOURCES}
)

# Include Directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
)

target_compile_options(${PROJECT_NAME} PRIVATE -w) # fix warning

# Link Libraries
target_link_libraries(${PROJECT_NAME}
    Boost::boost
    bdsg
    handlegraph
    Eigen3::Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
    gbwt
)

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
    src/snarl_analyser.cpp
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
