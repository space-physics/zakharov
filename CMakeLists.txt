cmake_minimum_required(VERSION 3.0)
project(zahkarov CXX Fortran)
enable_testing()

include(cmake/compilers.cmake)

add_executable(zakhfort zakh.f90 perf.f90)
target_compile_options(zakhfort PRIVATE ${FFLAGS})

#------------
find_package(Boost COMPONENTS filesystem program_options)
if(Boost_FOUND)
  add_executable(zakh zakh.cpp)
  target_link_libraries(zakh PRIVATE Boost::filesystem Boost::program_options)
endif()

