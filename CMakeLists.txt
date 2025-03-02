cmake_minimum_required(VERSION 3.10)

project(stoat VERSION 0.1.1 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${MAKE_BINARY_DIR}/bin") 

# Enable testing
enable_testing()

include(CTest)
include(lib/Catch/Catch.cmake)

# Include directories
include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

# Find Boost
find_package(Boost)

# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})
include_directories ("${PROJECT_SOURCE_DIR}/lib/Catch")
include_directories(BEFORE SYSTEM "${PROJECT_SOURCE_DIR}")

# Link libraries
target_link_libraries(${PROJECT_NAME} 
    Boost::boost 
    bdsg 
    handlegraph 
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
)

# Add the source directory
add_subdirectory(src)
add_subdirectory(tests/unittest)
add_subdirectory(tests/fonctionnal)
