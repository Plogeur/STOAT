project(stoat_fonctionnal)

# Find Boost
find_package(Boost)

# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Add test executable
set(TEST_SRCS
    binary_test.cpp
    quantitative_test.cpp
)

# Register test with CTest
add_test(NAME Stoat_test COMMAND run_tests_fonctionnal)

# Create test executable
add_executable(run_tests_fonctionnal binary_test.cpp quantitative_test.cpp)

# Link required libraries (if needed)
target_link_libraries(run_tests_fonctionnal PRIVATE)
