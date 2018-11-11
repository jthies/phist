  if(PHIST_HAVE_MKL)
    set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -DMKL_LP64")
    set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -DMKL_LP64")
  endif()

  if(NOT PHIST_CROSS_COMPILE)
    set(ARCH_FLAG "-march=native")
  endif()

  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # FLANG is still experimental and one may want to use gfortran instead
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  endif()

  # -ffast-math kills high precision stuff!
  set(FAST_MATH "-fno-math-errno -ffinite-math-only -fno-signed-zeros -fno-trapping-math")
  set (CMAKE_C_FLAGS_RELEASE        "-O3 ${FAST_MATH} ${ARCH_FLAG}")
  set (CMAKE_CXX_FLAGS_RELEASE      "-O3 ${FAST_MATH} ${ARCH_FLAG}")
  set (CMAKE_Fortran_FLAGS_RELEASE  "-O3 ${FAST_MATH} ${ARCH_FLAG}")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -Wall")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -Wall")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall")

  if (CMAKE_Fortran_COMPILER_ID STREQUAL "Clang")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_Debug} -Mbounds")
  endif()
