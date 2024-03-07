cmake_minimum_required(VERSION 3.27)
project(second C)

set(CMAKE_C_STANDARD 99)

add_executable(second
        main.c
        sequential_program.h
        sequential_program.c
        parallel_program.h
        parallel_program.c
)

target_link_libraries(second PRIVATE m)