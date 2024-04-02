cmake_minimum_required(VERSION 3.27)
project(first C)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 23)

set(SOURCES main.c
        parallel_program.c
        parallel_program.h
)

add_executable(first main.c
        parallel_program.h
        parallel_program.c
        program.c
        program.h)

target_include_directories(first PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(first PRIVATE ${MPI_LIBRARIES} m)
