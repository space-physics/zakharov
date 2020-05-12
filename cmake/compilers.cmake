if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native -Wall -Wextra -fexceptions -Werror=array-bounds)

  string(APPEND CMAKE_Fortran_FLAGS " -fimplicit-none")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fcheck=all -ffpe-trap=invalid,zero,overflow")

  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 8)
    string(APPEND CMAKE_Fortran_FLAGS " -std=f2018")
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  if(WIN32)
    string(APPEND CMAKE_Fortran_FLAGS " /warn /traceback")
  else()
    string(APPEND CMAKE_Fortran_FLAGS " -warn declarations -warn -traceback")
  endif()

  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fpe0 -debug extended -check all")
endif()
