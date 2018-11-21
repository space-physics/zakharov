cmake_minimum_required(VERSION 3.7)
project(zahkarov CXX Fortran)
enable_testing()

include(cmake/compilers.cmake)

add_executable(zakhfort zakh.f90 perf.f90)
target_compile_options(zakhfort PRIVATE ${FFLAGS})
target_link_libraries(zakhfort ${OpenCoarrays_LIBRARIES})

#------------
find_package(Boost COMPONENTS system filesystem program_options)
if(Boost_FOUND)
  add_executable(zakh zakh.cpp)
  target_link_libraries(zakh PRIVATE ${Boost_LIBRARIES})
  target_include_directories(zakh PRIVATE ${Boost_INCLUDE_DIRS})
endif()

