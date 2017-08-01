cmake_minimum_required (VERSION 2.8.12)
project(zahkarov CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3 -mtune=native -Wall -Wextra -Wpedantic -fexceptions -Warray-bounds -fstack-protector)

#------------
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIR})


add_executable(zakh zakh.cpp)
target_link_libraries(zakh -lboost_system -lboost_filesystem -lboost_program_options)

add_executable(zakhfort zakh.f90)
