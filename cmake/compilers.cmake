
set(CMAKE_CXX_STANDARD 17)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)

  add_compile_options(-mtune=native -Wall -Wextra -Wpedantic -fexceptions -Werror=array-bounds)

  list(APPEND FFLAGS -Warray-temporaries -Wconversion -fimplicit-none)

  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    list(APPEND FFLAGS -fcheck=all -ffpe-trap=invalid,zero,overflow)
  endif()

  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
    list(APPEND FFLAGS -std=f2018)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  if(WIN32)
    list(APPEND FFLAGS /4Yd /warn /traceback /coarray)
  else()
    list(APPEND FFLAGS -implicitnone -warn -traceback -coarray)
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    list(APPEND FFLAGS -fpe0 -debug extended -check all)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL PGI)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Flang)

endif()
