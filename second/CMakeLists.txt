cmake_minimum_required(VERSION 3.27)
project(second C)

set(CMAKE_C_STANDARD 99)

set(SOURCES
        main.c
        sequential_program.h
        sequential_program.c
        parallel_for_program.h
        parallel_for_program.c
        parallel_section_program.h
        parallel_section_program.c
        schedule_program.h
        schedule_program.c
)

add_executable(main.out ${SOURCES})

set(CMAKE_C_COMPILER gcc)

find_package(OpenMP REQUIRED)
if(OpenMP_C_FOUND)
    target_link_libraries(main.out PUBLIC OpenMP::OpenMP_C)
endif()

target_link_libraries(main.out PRIVATE m)
