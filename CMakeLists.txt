cmake_minimum_required(VERSION 3.15)
project(graph-algorithms VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

include_directories(include)

add_library(graph-core
    src/Graph.cpp
    src/HungarianAlgorithm.cpp
    src/GraphGenerator.cpp
)

target_include_directories(graph-core PUBLIC include)

add_executable(graph-generator
    src/main_generator.cpp
)

target_link_libraries(graph-generator graph-core)

add_executable(graph-solver
    src/main_solver.cpp
)

target_link_libraries(graph-solver graph-core)

install(TARGETS graph-generator graph-solver DESTINATION bin)
