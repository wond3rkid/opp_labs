cmake_minimum_required(VERSION 3.27)
project(secondparallel C)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 23)

set(SOURCES main.c)

add_executable(secondparallel main.c)

target_include_directories(secondparallel PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(secondparallel PRIVATE ${MPI_LIBRARIES} m)