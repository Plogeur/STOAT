cmake_minimum_required(VERSION 3.10)
project(stoat_cxx VERSION 0.1.1 LANGUAGES CXX)

# C++ Standard Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Find Packages
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# Find Python
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Ensure dependencies are installed
add_custom_target(install_python_dependencies ALL
    COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade pip
    COMMAND ${Python3_EXECUTABLE} -m pip install numpy pandas matplotlib scipy qmplot statsmodels
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Add Executable
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/utils.cpp
    src/arg_parser.cpp
    src/binary_analysis.cpp
    src/adjusted_pvalue.cpp
    src/gaf_creator.cpp
    src/list_snarl_paths.cpp
    src/matrix.cpp
    src/quantitative_analysis.cpp
    src/snarl_parser.cpp
    src/lmm.cpp
    ${ALGLIB_SOURCES}
)

add_dependencies(${PROJECT_NAME} install_python_dependencies)

# Include Directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${Boost_INCLUDE_DIRS}
    ${HTSLIB_INCLUDE_DIRS}
)

# Link Libraries
target_link_libraries(${PROJECT_NAME}
    Boost::boost
    bdsg
    handlegraph
    Eigen3::Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
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
    tests/unittest/binary_analysis_unit.cpp
    tests/unittest/gaf_creator_unit.cpp
    tests/unittest/list_snarl_paths_unit.cpp
    tests/unittest/matrix_unit.cpp
    tests/unittest/quantitative_analysis_unit.cpp

    src/utils.cpp
    src/arg_parser.cpp
    src/binary_analysis.cpp
    src/gaf_creator.cpp
    src/list_snarl_paths.cpp
    src/matrix.cpp
    src/quantitative_analysis.cpp
    src/snarl_parser.cpp
    src/lmm.cpp
)

# add src target lib + Catch2 lib 
target_link_libraries(unit_tests PRIVATE 
    Catch2::Catch2WithMain     
    Boost::boost
    bdsg
    handlegraph
    Eigen3::Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES})
    