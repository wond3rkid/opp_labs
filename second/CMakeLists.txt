cmake_minimum_required(VERSION 3.27)
project(second C)

set(CMAKE_C_STANDARD 11)

add_executable(second
        main.c
        process.h
        process.c
)
