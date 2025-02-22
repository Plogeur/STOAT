cmake_minimum_required(VERSION 3.10)

# Project name and version
project(stoat_cxx VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compile options (equivalent to CXXFLAGS)
add_compile_options(-Wall -Wextra -O2)

# Include directories
include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

# Find Boost
find_package(Boost REQUIRED)

# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Source files
set(SRCS 
    src/matrix.cpp 
    src/arg_parser.cpp 
    src/snarl_parser.cpp 
    src/binary_analysis.cpp 
    src/quantitative_analysis.cpp 
    src/list_snarl_paths.cpp
    src/main.cpp
)

# Create the executable
add_executable(${PROJECT_NAME} ${SRCS})

# Link libraries
target_link_libraries(${PROJECT_NAME} 
    Boost::boost 
    bdsg 
    handlegraph 
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
)

# Ensure proper linking flags for HTSlib
target_compile_options(${PROJECT_NAME} PRIVATE ${HTSLIB_CFLAGS_OTHER})
