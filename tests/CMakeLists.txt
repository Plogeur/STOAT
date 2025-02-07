cmake_minimum_required(VERSION 3.10)
project(list_snarl_paths)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/include /usr/local/include)
link_directories(/usr/include /usr/local/lib)  # libbdsg.dylib

add_executable(list_snarl_paths src/main.cpp src/list_snarl_paths.cpp)

# Link the bdsg library to your executable target
target_link_libraries(list_snarl_paths bdsg handlegraph libsdsl.so.3)
