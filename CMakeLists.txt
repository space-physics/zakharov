cmake_minimum_required(VERSION 3.12)

project(zahkarov
  LANGUAGES CXX Fortran)

include(CTest)

include(cmake/compilers.cmake)

add_executable(zakhfort zakh.f90 perf.f90)
add_test(NAME ZakharovFortran COMMAND $<TARGET_FILE:zakhfort> ${CMAKE_CURRENT_BINARY_DIR} 1e-6 50 100 250)

find_package(Boost COMPONENTS filesystem program_options)
if(Boost_FOUND)
  add_executable(zakh zakh.cpp)
  target_link_libraries(zakh PRIVATE Boost::filesystem Boost::program_options)
endif()
