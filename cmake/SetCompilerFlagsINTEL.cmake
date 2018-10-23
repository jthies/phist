  set (CMAKE_C_FLAGS        "-align -Wno-unused-variable ${MKL_FLAG}")
  set (CMAKE_CXX_FLAGS      "-align -Wno-unused-variable ${MKL_FLAG}")
  set (CMAKE_Fortran_FLAGS  "${MKL_FLAG} -align array64byte")

  set (CMAKE_C_FLAGS_RELEASE        "-O2 -prec-div -xHOST")
  set (CMAKE_CXX_FLAGS_RELEASE      "-O2 -prec-div -xHOST")
  set (CMAKE_Fortran_FLAGS_RELEASE  "-O3 -prec-div -xHOST")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -debug")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -debug")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -debug")

#  set (CMAKE_C_FLAGS_DEBUG        "-check-pointers=rw -fstack-protector -w2 -debug -traceback")
#  set (CMAKE_CXX_FLAGS_DEBUG      "-check-pointers=rw -w2 -fstack-protector -debug -traceback")
  set (CMAKE_C_FLAGS_DEBUG        "-fstack-protector -w2 -debug -traceback")
  set (CMAKE_CXX_FLAGS_DEBUG      "-fstack-protector -w2 -debug -traceback")
  set (CMAKE_Fortran_FLAGS_DEBUG  "-check all -check noarg_temp_created -debug -traceback")

  # this is for some reason required with Intel 15 on Emmy:
  set (CMAKE_EXE_LINKER_FLAGS_DEBUG "-lifcore -lifport")

  # for ipo
  # melven: probably not useful or even dangerous?
  #  kernel routines should probably best optimized individually;
  #  for builtin kernels: I'm passing arrays as scalars in builtin kernels (call by reference),
  #               because the kernel routines would be too complex otherwise
  #               (and they should be fast and not complex!).
  #               The compiler best shouldn't be able to know about this...
  #set_property(DIRECTORY . PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
