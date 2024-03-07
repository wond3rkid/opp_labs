cmake_minimum_required(VERSION 3.27)
project(second C)

set(CMAKE_C_STANDARD 99)
find_package(OpenMP REQUIRED)

# Добавление компиляторных флагов OpenMP
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(second
        main.c
        sequential_program.h
        sequential_program.c
        parallel_program.h
        parallel_program.c
)

target_link_libraries(second PRIVATE m)
target_link_libraries(second PUBLIC OpenMP::OpenMP_C)