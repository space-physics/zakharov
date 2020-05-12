cmake_minimum_required(VERSION 3.12)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(zahkarov
  LANGUAGES CXX Fortran)
enable_testing()

include(cmake/compilers.cmake)

add_executable(zakhfort zakh.f90 perf.f90)
add_test(NAME ZakharovFortran COMMAND $<TARGET_FILE:zakhfort> ${CMAKE_CURRENT_BINARY_DIR} 1e-6 50 100 250)

find_package(Boost COMPONENTS filesystem program_options)
if(Boost_FOUND)
  add_executable(zakh zakh.cpp)
  target_link_libraries(zakh PRIVATE Boost::filesystem Boost::program_options)
endif()
