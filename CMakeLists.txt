cmake_minimum_required (VERSION 2.8.12)
project(zahkarov CXX)

add_compile_options(-mtune=native -Wall -Wextra -Wpedantic -fexceptions -Warray-bounds)

add_executable(zakh zakh.cpp)
set_target_properties(zakh PROPERTIES CXX_STANDARD 14)
#------------
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(zakh ${Boost_LIBRARIES})

