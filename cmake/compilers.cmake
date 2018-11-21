if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

set(CMAKE_CXX_STANDARD 17)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

  find_package(OpenCoarrays)
  if(OpenCoarrays_FOUND)
    list(APPEND OpenCoarrays_LIBRARIES caf_mpi)
    list(APPEND FFLAGS -fcoarray=lib)
  else()
    list(APPEND FFLAGS -fcoarray=single)
  endif()

  add_compile_options(-mtune=native -Wall -Wextra -Wpedantic -fexceptions -Werror=array-bounds)

  list(APPEND FFLAGS -Warray-temporaries -Wconversion -fimplicit-none 
-fcheck=all -ffpe-trap=invalid,zero,overflow)


elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
  list(APPEND FFLAGS -warn -fpe0 -traceback -debug extended)# -check all -coarray )
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL Flang)  
  list(APPEND FFLAGS -Mallocatable=03)
  list(APPEND FLIBS -static-flang-libs)
endif()
