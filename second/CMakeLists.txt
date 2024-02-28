cmake_minimum_required(VERSION 3.27)
project(second C)

set(CMAKE_C_STANDARD 99)

add_executable(second
        main.c
        process.h
        process.c
)
#m - math.h
target_link_libraries(second m)