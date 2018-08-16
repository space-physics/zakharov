cmake_minimum_required (VERSION 3.7)
project(zahkarov CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

add_compile_options(-O3 -g)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(zakhfort)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

  find_package(OpenCoarrays)
  if(OpenCoarrays_FOUND)
    list(APPEND OpenCoarrays_LIBRARIES caf_mpi)
    add_compile_options(-fcoarray=lib)
  else()
    add_compile_options(-fcoarray=single)
  endif()

  add_compile_options(-mtune=native -Wall -Wextra -Wpedantic -fexceptions -Werror=array-bounds)

  target_compile_options(zakhfort PUBLIC 
  -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace 
  -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow)

  target_link_libraries(zakhfort ${OpenCoarrays_LIBRARIES})

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  target_compile_options(zakhfort PUBLIC -warn -fpe0 -traceback -debug extended)# -check all -coarray )
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)  
  target_compile_options(zakhfort PUBLIC -Mallocatable=03)
  list(APPEND FFLAGS -static-flang-libs)
endif()

#------------
find_package(Boost COMPONENTS system filesystem program_options)
if(Boost_FOUND)
  add_executable(zakh zakh.cpp)
  target_link_libraries(zakh PRIVATE ${Boost_LIBRARIES})
  target_include_directories(zakh PRIVATE ${Boost_INCLUDE_DIRS})
endif()
#------------
target_sources(zakhfort PRIVATE zakh.f90 perf.f90)
target_link_libraries(zakhfort PRIVATE  ${FFLAGS})

