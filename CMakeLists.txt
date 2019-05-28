cmake_minimum_required(VERSION 3.0)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()
project(zahkarov
  LANGUAGES NONE)
enable_testing()

include(cmake/compilers.cmake)

include(CheckLanguage)
check_language(Fortran)

if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)

  add_executable(zakhfort zakh.f90 perf.f90)
  target_compile_options(zakhfort PRIVATE ${FFLAGS})
endif()


check_language(CXX)
if(CMAKE_CXX_COMPILER)
  enable_language(CXX)

  set(CMAKE_CXX_STANDARD 17)

  find_package(Boost COMPONENTS filesystem program_options)
  if(Boost_FOUND)
    add_executable(zakh zakh.cpp)
    target_link_libraries(zakh PRIVATE Boost::filesystem Boost::program_options)
  endif()
endif()
