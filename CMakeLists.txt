cmake_minimum_required(VERSION 3.10)
project(stoat_cxx VERSION 0.1.1 LANGUAGES CXX)

# C++ Standard Settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Optional: Set runtime output directory
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Find Packages
find_package(Boost REQUIRED)
find_package(Eigen REQUIRED)

# Add Executable
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/arg_parser.cpp
    src/binary_analysis.cpp
    src/gaf_creator.cpp
    src/list_snarl_paths.cpp
    src/matrix.cpp
    src/quantitative_analysis.cpp
    src/snarl_parser.cpp
    ${ALGLIB_SOURCES}
)

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
    Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
)

# Enable testing
# enable_testing()
# add_subdirectory(tests/unittest)

# include(CTest)
# include(lib/Catch/Catch.cmake)

# # Include directories
# include_directories(/usr/include /usr/local/include)
# link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

# # Find HTSlib
# find_package(PkgConfig REQUIRED)
# pkg_check_modules(HTSLIB REQUIRED htslib)

# include_directories(${HTSLIB_INCLUDE_DIRS})
# link_directories(${HTSLIB_LIBRARY_DIRS})

# include_directories ("${PROJECT_SOURCE_DIR}/lib/Catch")
# include_directories(BEFORE SYSTEM "${PROJECT_SOURCE_DIR}")

# # Add the executable
# add_executable(${PROJECT_NAME} 
#     src/main.cpp
#     src/arg_parser.cpp
#     src/binary_analysis.cpp
#     src/gaf_creator.cpp
#     src/list_snarl_paths.cpp
#     src/matrix.cpp
#     src/quantitative_analysis.cpp
#     src/snarl_parser.cpp             
# )

# # Link libraries
# target_link_libraries(${PROJECT_NAME} 
#     Boost::boost 
#     bdsg 
#     handlegraph 
#     libsdsl.so.3
#     ${HTSLIB_LIBRARIES}
# )

# # Add the test executable
# add_executable(TestRunner 
#     tests/unittest/binary_analysis_unit.cpp
#     # tests/unittest/gaf_creator_unit.cpp
#     # tests/unittest/list_snarl_paths_unit.cpp
#     tests/unittest/matrix_unit.cpp
#     tests/unittest/quantitative_analysis_unit.cpp

#     src/arg_parser.cpp
#     src/binary_analysis.cpp
#     src/gaf_creator.cpp
#     src/list_snarl_paths.cpp
#     src/matrix.cpp
#     src/quantitative_analysis.cpp
#     src/snarl_parser.cpp   
# )

# # Ensure it finds the headers from src/
# target_include_directories(TestRunner PRIVATE "${CMAKE_SOURCE_DIR}/src")

# # Link necessary libraries
# target_link_libraries(TestRunner PRIVATE Boost::boost bdsg handlegraph libsdsl.so.3 Catch2 ${HTSLIB_LIBRARIES})

# # Register the test with CTest
# add_test(NAME MyTests COMMAND TestRunner)
