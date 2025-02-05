cmake_minimum_required(VERSION 3.10)
project(Snarl)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include)  # bdsg is in /usr/local/include
link_directories(/usr/local/lib)  # libbdsg.dylib is in /usr/local/lib

add_executable(Snarl src/main.cpp src/list_snarl_paths.cpp)

# Link the bdsg library to your executable target
target_link_libraries(Snarl bdsg)