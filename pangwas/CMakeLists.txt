
# CMake version
cmake_minimum_required(VERSION 3.15)

# Project name
project(pangwas)

include(GNUInstallDirs)

# c++14
set(CMAKE_CXX_STANDARD 14)


option(OPTIMIZE "Build with optimization" ON)

if (OPTIMIZE)
# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif ()


# Find our pangwas package directory where input sources and dependencies are
set(src_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(deps_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps")

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)


# libbdsg (full build using its cmake config)
# We need it for a few overlays and overlay helpers used for the IntegratedSnarlFinder.
# Produces libbdsg target, dynamic only.
# Comes along with its *own* libhandlegraph, which it depends on.
# Make sure to never build its Python bindings.
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings for libbdsg" FORCE)
set(RUN_DOXYGEN OFF CACHE BOOL "Run libbdsg Doxygen" FORCE)
add_subdirectory(${deps_DIR}/libbdsg)


# Add libvgio using its cmake build
add_subdirectory(${deps_DIR}/libvgio)

# Copied from https://github.com/vgteam/vg/wiki/Programming-with-the-vg-API
# SDSL library ("sdsl" target) already exists from libbdsg

# GBWT and GBWTGraph don't have CMake yet so we need to explain how to build them.

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
find_package( OpenMP REQUIRED )

# GBWT
foreach (gbwt_source IN ITEMS algorithms.cpp bwtmerge.cpp cached_gbwt.cpp dynamic_gbwt.cpp fast_locate.cpp files.cpp gbwt.cpp internal.cpp metadata.cpp support.cpp test.cpp utils.cpp variants.cpp)
    list(APPEND gbwt_libFiles "deps/gbwt/src/${gbwt_source}")
endforeach()
if(MSVC)
    file(GLOB gbwt_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${deps_DIR}/gbwt/include/gbwt/*.h")
endif()
set( gbwt_SRCS ${gbwt_libFiles} ${gbwt_headerFiles})
add_library( gbwt ${gbwt_SRCS} )
target_include_directories( gbwt PUBLIC "${deps_DIR}/gbwt/include")
target_link_libraries( gbwt PUBLIC sdsl OpenMP::OpenMP_CXX )
install(TARGETS gbwt
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
 
# GBWTGraph
foreach (gbwtgraph_source IN ITEMS algorithms.cpp cached_gbwtgraph.cpp gbwtgraph.cpp gbz.cpp gfa.cpp internal.cpp minimizer.cpp path_cover.cpp utils.cpp)
    list(APPEND gbwtgraph_libFiles "deps/gbwtgraph/src/${gbwtgraph_source}")
endforeach()
if(MSVC)
    file(GLOB gbwtgraph_headerFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${deps_DIR}/gbwtgraph/include/gbwt/*.h")
endif()
set( gbwtgraph_SRCS ${gbwtgraph_libFiles} ${gbwtgraph_headerFiles})
add_library( gbwtgraph ${gbwtgraph_SRCS} )
target_include_directories( gbwtgraph PUBLIC "${deps_DIR}/gbwtgraph/include")
target_link_libraries( gbwtgraph PUBLIC gbwt handlegraph_shared )
install(TARGETS gbwtgraph
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
        
add_library(pangwas_objs OBJECT
    src/utils.cpp
    src/association_finder.cpp
    src/path_association_finder.cpp
    src/tester.cpp
    src/io/register_loader_saver_gbz.cpp
    src/io/register_loader_saver_gbzgraph.cpp
    src/io/register_loader_saver_hash_graph.cpp
    src/io/register_loader_saver_packed_graph.cpp
    src/io/register_io.cpp
)

target_include_directories(pangwas_objs PUBLIC Boost::math)

target_link_libraries(pangwas_objs libbdsg
                                   gbwtgraph
                                   VGio::vgio)

add_executable(pangwas ${src_DIR}/main.cpp)

target_link_libraries(pangwas pangwas_objs)


# Make tests
add_executable(test_pangwas
    src/test/main_test.cpp
    src/test/test_association_finder.cpp
    src/test/test_path_association_finder.cpp
    src/test/test_tester.cpp
)
target_link_libraries(test_pangwas pangwas_objs)
# Add catch2 for unittesting
include_directories(${deps_DIR}/catch2/single_include/catch2 )
