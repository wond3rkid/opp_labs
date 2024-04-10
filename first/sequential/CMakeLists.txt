cmake_minimum_required(VERSION 3.27)
project(sequential C)

set(CMAKE_C_STANDARD 11)

add_executable(sequential
        main.c)

target_link_libraries(sequential  m)