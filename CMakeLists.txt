cmake_minimum_required (VERSION 2.8.12)
project(zahkarov CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-mtune=native -Wall -Wextra -Wpedantic -fexceptions -Warray-bounds)

#------------
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIRS})


add_executable(zakh zakh.cpp)
target_link_libraries(zakh ${Boost_LIBRARIES})

