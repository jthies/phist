string(TOUPPER ${PHIST_KERNEL_LIB} PHIST_KERNEL_LIB_CAPS)
set(PHIST_KERNEL_LIB_${PHIST_KERNEL_LIB_CAPS} ON)

if  (PHIST_KERNEL_LIB STREQUAL "epetra")

elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")

elseif (PHIST_KERNEL_LIB STREQUAL "ghost")
  # we already link with ghost because of the queuing system  
  set(PHIST_USE_GHOST ON)
  set(PHIST_MVECS_ROW_MAJOR 1 CACHE STRING "wether multi-vectors should be stored in row-major ordering")
  set(PHIST_USE_GHOST_TASKS 0 CACHE STRING "execute all kernels as separate GHOST tasks, required for using PHISTs tasking features.")

  set(PHIST_SELL_C -1 CACHE STRING "chunk size C for the SELL-C-sigma format (-1: automatic choice)")
  set(PHIST_SELL_SIGMA -1 CACHE STRING "sorting scope sigma for the SELL-C-sigma format (-1: automatic choice)")

  set(PHIST_DEFAULT_SPMV_MODE 0 CACHE STRING "default GHOST spmv mode (0: VECTOR, 8: OVERLAP, 16: TASK)")

elseif (PHIST_KERNEL_LIB STREQUAL "builtin")
  set(PHIST_MVECS_ROW_MAJOR 1)
  option(PHIST_BUILTIN_MEMALIGN "Specify desired alignment for builtin Fortran kernels (in byte). Setting it to 0 will use regular mallocs." 0)
  if(PHIST_HAVE_AVX2)
    option(PHIST_HIGH_PRECISION_KERNELS "Activate experimental high-precision builtin kernel routines" Off)
    option(PHIST_HIGH_PRECISION_KERNELS_FORCE "Always force the usage of high precision reductions even if this requires copying data and is slow!" Off)
  endif()

  if( NOT PHIST_ENABLE_MPI )
    message(FATAL_ERROR "Builtin kernel lib requires MPI!")
  endif()
elseif (PHIST_KERNEL_LIB STREQUAL "magma")
  set(PHIST_MVECS_ROW_MAJOR 0)
elseif (PHIST_KERNEL_LIB STREQUAL "petsc")
  set(PHIST_MVECS_ROW_MAJOR 0)
elseif (PHIST_KERNEL_LIB STREQUAL "eigen")
  set(PHIST_MVECS_ROW_MAJOR 0 CACHE STRING "wether multi-vectors should be stored in row-major ordering")
  if( PHIST_HAVE_MPI )
    message(FATAL_ERROR "The eigen kernel lib has no MPI support!")
  endif()
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set or not recognized" )
endif ()

