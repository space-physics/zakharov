cmake_minimum_required (VERSION 2.8.12)
project(zahkarov CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3 -mtune=native -g -Wall -Wextra -Wpedantic -fexceptions -Warray-bounds -fstack-protector)

#------------
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
include_directories(${Boost_INCLUDE_DIR})


add_executable(zakh zakh.cpp)
target_link_libraries(zakh -lboost_system -lboost_filesystem -lboost_program_options)

add_executable(zakhfort zakh.f90)
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    target_compile_options(zakhfort PUBLIC 
    -std=f2008
    -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
    target_compile_options(zakhfort PUBLIC -check all -fpe0 -warn -traceback -debug extended)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")  # FLANG
    target_compile_options(zakhfort PUBLIC -Mallocatable=03)
endif()
