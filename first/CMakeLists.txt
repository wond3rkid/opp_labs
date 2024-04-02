cmake_minimum_required(VERSION 3.27)
project(first C)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 23)

set(SOURCES main.c)

add_executable(first main.c)

target_include_directories(first PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(first PRIVATE ${MPI_LIBRARIES} m)
