cmake_minimum_required(VERSION 3.10)
project(slink VERSION 0.1.0 LANGUAGES CXX)

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

# Add Executable
add_executable(${PROJECT_NAME}
    src/main.cpp
    src/arg_parser.cpp
    src/matrix.cpp
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
    Eigen3::Eigen
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
)
