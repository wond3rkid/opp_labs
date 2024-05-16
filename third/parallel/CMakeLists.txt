cmake_minimum_required(VERSION 3.27)
project(third C)

find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 23)

set(SOURCES main.c)

add_executable(third main.c)

target_include_directories(third PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(third PRIVATE ${MPI_LIBRARIES} m)