cmake_minimum_required (VERSION 2.8.12)
project(zahkarov CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_compile_options(-O3 -mtune=native -g)

#------------
if(NOT WIN32)
  find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
  include_directories(${Boost_INCLUDE_DIR})

  add_executable(zakh zakh.cpp)
  target_link_libraries(zakh ${Boost_LIBRARIES})
endif()
#------------
add_executable(zakhfort zakh.f90 perf.f90)
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

  find_package(CAF)
  if(CAF_FOUND)
    set(FCOOPT -fcoarray=lib)
    set(FCOLIB caf_mpi)
  else()
    set(FCOOPT -fcoarray=single)
  endif()

  add_compile_options(-Wall -Wextra -Wpedantic -fexceptions -Werror=array-bounds)

  target_compile_options(zakhfort PUBLIC 
  ${FCOOPT}
  -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace 
  -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow)

  target_link_libraries(zakhfort ${FCOLIB})


elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  target_compile_options(zakhfort PUBLIC -warn -fpe0 -traceback -debug extended)# -check all -coarray )
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")  # FLANG
  target_compile_options(zakhfort PUBLIC -Mallocatable=03)
endif()
