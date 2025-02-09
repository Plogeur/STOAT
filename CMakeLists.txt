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

# Define Boost dependency
find_package(Boost REQUIRED COMPONENTS math)

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Boost not found. Please install Boost and try again.")
endif()

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

# Link Boost libraries
target_link_libraries(${PROJECT_NAME} Boost::math bdsg handlegraph libsdsl.so.3)

# Custom target for get-deps
add_custom_target(get-deps
    COMMAND ${CMAKE_COMMAND} -E echo "Fetching dependencies..."
    COMMAND ${CMAKE_COMMAND} -E echo "Boost is required and should be pre-installed."
    COMMAND ${CMAKE_COMMAND} -E echo "On Debian-based systems, run: sudo apt install libboost-all-dev"
    COMMAND ${CMAKE_COMMAND} -E echo "On macOS, run: brew install boost"
    COMMAND ${CMAKE_COMMAND} -E echo "For source installation, see: https://www.boost.org"
)
