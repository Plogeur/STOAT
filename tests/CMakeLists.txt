cmake_minimum_required(VERSION 3.10)
project(Snarl)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

add_executable(Snarl src/main.cpp src/list_snarl_paths.cpp)

# Link the bdsg library to your executable target
target_link_libraries(Snarl bdsg)
