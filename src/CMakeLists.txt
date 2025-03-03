project(stoat_src)

aux_source_directory(. SOURCES_LIST)
# static lib are include in the final bin
add_library(${PROJECT_NAME} STATIC ${SOURCES_LIST})

set(SRC
    main.cpp
    arg_parser.cpp
    binary_analysis.cpp
    gaf_creator.cpp
    list_snarl_paths.cpp
    matrix.cpp
    quantitative_analysis.cpp
    snarl_parser.cpp
)

# Create main executable
add_executable(${PROJECT_NAME} 
    main.cpp
    arg_parser.cpp
    binary_analysis.cpp
    gaf_creator.cpp
    list_snarl_paths.cpp
    matrix.cpp
    quantitative_analysis.cpp
    snarl_parser.cpp    
)
