cmake_minimum_required (VERSION 2.8.12)
project(zahkarov CXX)

add_compile_options(-std=c++11 -mtune=native -ffast-math -Wall -Wextra -Wpedantic -fexceptions -Warray-bounds)

add_executable(zakh zakh.cpp)
#------------
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(zakh ${Boost_LIBRARIES})

