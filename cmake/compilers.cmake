if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-mtune=native -Wall -Wextra -fexceptions -Werror=array-bounds)

  string(APPEND CMAKE_Fortran_FLAGS " -fimplicit-none")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fcheck=all -ffpe-trap=invalid,zero,overflow")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

  if(WIN32)
    string(APPEND CMAKE_Fortran_FLAGS " /warn /traceback")
  else()
    string(APPEND CMAKE_Fortran_FLAGS " -warn -traceback")
  endif()

  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fpe0 -debug extended -check all")
endif()
