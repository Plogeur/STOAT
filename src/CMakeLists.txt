# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)

include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

target_link_libraries(stoat_exe
    Boost::boost 
    bdsg 
    handlegraph 
    libsdsl.so.3
    ${HTSLIB_LIBRARIES}
)

add_executable(stoat_exe
    main.cpp
    arg_parser.cpp
    binary_analysis.cpp
    gaf_creator.cpp
    list_snarl_paths.cpp
    matrix.cpp
    quantitative_analysis.cpp
    snarl_parser.cpp             
)

target_include_directories(stoat_exe PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
