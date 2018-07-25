  set (CMAKE_C_FLAGS        "-std=c99")
  set (CMAKE_CXX_FLAGS      "-std=c++11")
  # assume we are using the LLVM Flang compiler (PGI flavour)
  set (CMAKE_Fortran_FLAGS  "-Mpreprocess -DFLANG")

  if(PHIST_HAVE_MKL)
    set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -DMKL_LP64")
    set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -DMKL_LP64")
  endif()

  # -ffast-math kills high precision stuff!
  set(FAST_MATH "-fno-math-errno -ffinite-math-only -fno-signed-zeros -fno-trapping-math")
  set (CMAKE_C_FLAGS_RELEASE        "-O3 ${FAST_MATH} -march=native")
  set (CMAKE_CXX_FLAGS_RELEASE      "-O2")
  set (CMAKE_Fortran_FLAGS_RELEASE  "-march=native -O3 ${FAST_MATH}")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -Wall")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -Wall")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Mbounds")
