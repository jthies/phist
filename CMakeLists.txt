cmake_minimum_required(VERSION 2.8.10)

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Trilinos")

# Define some custom variables
project(phist C CXX Fortran)
SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)
set(PHIST_VERSION 0.8.2)
set(PRODUCT_NAME phist-${PHIST_VERSION})

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Release, RelWithDebug, Debug" FORCE)
else()
  set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: Release, RelWithDebug, Debug" FORCE)
endif()

if (NOT BUILD_SHARED_LIBS)
    message(STATUS "No library type specified. Setting it to shared.")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Set to ON/OFF for shared/static library build.")
endif()


option(PHIST_ENABLE_MPI "Enable MPI within PHIST" ON)
option(PHIST_ENABLE_OPENMP "Enable OpenMP within PHIST" ON)
option(PHIST_ENABLE_SP "Enable single precision functions and drivers" OFF)
option(PHIST_ENABLE_COMPLEX_TESTS "Enable unit tests for complex data types" ON)

########################################################################
# Lists of libraries, executables, tests
########################################################################
set(library_list)
set(driver_list)
set(test_list)

########################################################################
# List of external libraries for executables and tests to link, e.g.
# GHOST, LAPACK. TODO: also for essex-physics?
########################################################################
set(external_libs)
set(PHIST_PKG_CONFIG_REQUIRE)
set(PHIST_PKG_CONFIG_CFLAGS)

########################################################################
# Deal with LAPACK and LAPACKe
# Specify paths in CMAKE_PREFIX_PATH.
########################################################################
find_path(LAPACKE_INCL NAMES "lapacke.h" "mkl_lapacke.h")
if(LAPACKE_INCL STREQUAL "LAPACKE_INCL-NOTFOUND")
  message(FATAL_ERROR "LAPACKe include file not found. Add the proper path to CMAKE_PREFIX_PATH.")
else()
  message(STATUS "LAPACKe include file found in: ${LAPACKE_INCL}")
  include_directories(${LAPACKE_INCL})
endif()

if(NOT CMAKE_C_COMPILER_ID STREQUAL "Intel")
# TODO: this assumes that if you use the intel compilers, you also have MKL (!)

  # LAPACKe
  unset(LAPACKE_LIBS CACHE)
  find_library(LAPACKE_LIBS NAMES "liblapacke.so" "liblapacke.a")
  if(LAPACKE_LIBS STREQUAL "LAPACKE_LIBS-NOTFOUND")
    message(FATAL_ERROR "LAPACKe libraries not found. Add the proper path to CMAKE_PREFIX_PATH.")
  else()
    message(STATUS "LAPACKe libraries found: ${LAPACKE_LIBS}")
    LIST(APPEND external_libs ${LAPACKE_LIBS})
  endif()

  # LAPACK
  unset(LAPACK_LIBS CACHE)
  find_library(LAPACK_LIBS NAMES "liblapack.so" "liblapack.a")
  if(LAPACK_LIBS STREQUAL "LAPACK_LIBS-NOTFOUND")
    message(FATAL_ERROR "LAPACK libraries not found. Add the proper path to CMAKE_PREFIX_PATH.")
  else()
    message(STATUS "LAPACK libraries found: ${LAPACK_LIBS}")
    LIST(APPEND external_libs ${LAPACK_LIBS})
  endif()

  # BLAS
  unset(BLAS_LIBS CACHE)
  find_library(BLAS_LIBS NAMES "libblas.so" "libblas.a")
  if(BLAS_LIBS STREQUAL "BLAS_LIBS-NOTFOUND")
    message(FATAL_ERROR "BLAS libraries not found. Add the proper path to CMAKE_PREFIX_PATH.")
  else()
    message(STATUS "BLAS libraries found: ${BLAS_LIBS}")
    LIST(APPEND external_libs ${BLAS_LIBS})
  endif()


endif()
########################################################################

# try to find a MKL header in the default location
if (NOT MKL_INCLUDE_DIR)
    message(STATUS "Try to find MKL in default paths and via MKLROOT")
    find_path(MKL_INCLUDE_DIR mkl.h $ENV{MKLROOT}/include)
endif()
if (MKL_INCLUDE_DIR)
    set(PHIST_HAVE_MKL 1)
    include_directories(${MKL_INCLUDE_DIR})
endif()


# output behavior
if(CMAKE_BUILD_TYPE MATCHES "Release")
	set(PHIST_OUTLEV 2 CACHE STRING "PHIST verbosity (release 2, 'INFO')")
else()
	set(PHIST_OUTLEV 3 CACHE STRING "PHIST verbosity (debug 3, 'VERBOSE')")
endif()

##############################################################################
# compiler-dependent flags 
##############################################################################


if (CMAKE_C_COMPILER_ID STREQUAL "GNU")

  message(STATUS "using GNU compiler flags")

  set (CMAKE_C_FLAGS        "-std=c99")
  set (CMAKE_CXX_FLAGS      "-std=c++0x")
  set (CMAKE_Fortran_FLAGS  "-cpp -ffree-line-length-none")
  if( PHIST_ENABLE_OPENMP )
    set(PHIST_HAVE_OPENMP 1)
    if (CMAKE_C_COMPILER_VERSION VERSION_LESS 4.9)
      set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -fopenmp")
      set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -fopenmp")
      set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fopenmp")
    else ()
      set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -fopenmp -fopenmp-simd")
      set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -fopenmp -fopenmp-simd")
      set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fopenmp -fopenmp-simd")
    endif ()
  endif()

  set (CMAKE_C_FLAGS_RELEASE        "-O3 -ffast-math -march=native")
  set (CMAKE_CXX_FLAGS_RELEASE      "-O2")
  set (CMAKE_Fortran_FLAGS_RELEASE  "-march=native -O3 -ffast-math")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")

  if (CMAKE_C_COMPILER_VERSION VERSION_LESS 4.8)
    set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter")
    set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter")
    set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -DTESTING -Wall -fcheck=all -fstack-protector-all")
  else ()
    set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -fstack-protector-all -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-local-typedefs")
    set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -fstack-protector-all -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-local-typedefs")
    set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -DTESTING -Wall -fcheck=all -fstack-protector-all")
  endif ()
  if (NOT CMAKE_C_COMPILER_VERSION VERSION_LESS 4.8)
    set(GCC_SANITIZE "" CACHE STRING "What should be sanitized by GCC for a Debug build? (address, thread, undefined, leak or leave empty)")
    if (GCC_SANITIZE)
      message(STATUS "Enabled ${GCC_SANITIZE}-sanitizer")
      if(GCC_SANITIZE STREQUAL "thread")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=${GCC_SANITIZE} -fPIC -pie")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${GCC_SANITIZE} -fPIC -pie")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=${GCC_SANITIZE} -fPIC -pie")
      else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=${GCC_SANITIZE}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=${GCC_SANITIZE}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=${GCC_SANITIZE}")
      endif()
    endif()
  endif()


elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")

  message(STATUS "using INTEL compiler flags")

  set (CMAKE_C_FLAGS        "-align -Wno-unused-variable -mkl -std=c99")
  set (CMAKE_CXX_FLAGS      "-align -Wno-unused-variable -mkl -std=c++0x")
  set (CMAKE_Fortran_FLAGS  "-fpp -mkl -align array64byte")

  if( PHIST_ENABLE_OPENMP )
    set(PHIST_HAVE_OPENMP 1)
    set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -openmp")
    set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -openmp")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -openmp")
  endif()

  set (CMAKE_C_FLAGS_RELEASE        "-O3 -xHOST -no-prec-div")
  set (CMAKE_CXX_FLAGS_RELEASE      "-O3 -xHOST -no-prec-div")
  set (CMAKE_Fortran_FLAGS_RELEASE  "-O3 -xHOST -no-prec-div")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -debug")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -debug")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -debug")

  set (CMAKE_C_FLAGS_DEBUG        "-check-pointers=rw -check=uninit,stack -fstack-protector -w2 -debug -traceback -DTESTING")
  set (CMAKE_CXX_FLAGS_DEBUG      "-check-pointers=rw -check=uninit,stack -fstack-protector -w2 -debug -traceback -DTESTING")
  set (CMAKE_Fortran_FLAGS_DEBUG  "-check all -debug -traceback -DTESTING")

  # for ipo
  # melven: probably not useful or even dangerous?
  #  kernel routines should probably best optimized individually;
  #  for builtin kernels: I'm passing arrays as scalars in builtin kernels (call by reference),
  #               because the kernel routines would be too complex otherwise
  #               (and they should be fast and not complex!).
  #               The compiler best shouldn't be able to know about this...
  #set_property(DIRECTORY . PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

else()

  message ("compiler ID: " ${CMAKE_C_COMPILER_ID})
  message ("No custom compiler flags are known.")

endif ()

# try to detect if we can use C++11 Lambdas
include(CheckCXXSourceCompiles)
include(CheckFortranSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
  int main(int argc, char** argv)
  {
    auto lfunc = [&](){ int _argc = argc; };
    return 0;
  }
" PHIST_HAVE_CXX11_LAMBDAS
)
if(PHIST_HAVE_CXX11_LAMBDAS)
  set(PHIST_PKG_CONFIG_CFLAGS "${PHIST_PKG_CONFIG_CFLAGS} -std=c++11")
endif()
CHECK_Fortran_SOURCE_COMPILES("
  program test
    implicit none
    integer :: i
    real :: x(8)
    !$omp simd
    do i = 1, 8
      x(i) = 0.
    end do
  end program test
" PHIST_HAVE_OPENMP_SIMD
)

##############################################################################
# set enable precision?
##############################################################################
if (${PHIST_ENABLE_SP})
  set(PHIST_HAVE_SP 1)
endif ()
##############################################################################
# setup MPI
##############################################################################

if (${PHIST_ENABLE_MPI})
  set(PHIST_HAVE_MPI 1)
#TODO - what is the best way to do this portably?
  find_package(MPI REQUIRED)
  if (${MPI_FOUND})
    set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    set (CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
    set (CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
    set (CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
    set (CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
# For now, just include everything
    include_directories (${MPI_INCLUDE_PATH})
    include_directories (${MPI_C_INCLUDE_PATH})
    include_directories (${MPI_CXX_INCLUDE_PATH})
    include_directories (${MPI_Fortran_INCLUDE_PATH})
  else()
    message(WARNING "could not find MPI. Trying to compile anyway, presuming the compiler/linker knows where to find it.")
  endif()
  set(MPIEXEC "mpirun" CACHE STRING "")
  set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")
endif()

##############################################################################
# setup which kernel library is being used.
##############################################################################

set(PHIST_KERNEL_LIB "builtin" CACHE STRING "Select underlying kernel libraries for phist (builtin|ghost|epetra|tpetra)")
message(STATUS "kernel lib " ${PHIST_KERNEL_LIB} " ...")

set(LOAD_TRILI_PKG,"")
if  (PHIST_KERNEL_LIB STREQUAL "epetra")
  list(APPEND LOAD_TRILI_PKG "EpetraExt")
  set(PHIST_KERNEL_LIB_EPETRA 1)
elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")
  list(APPEND LOAD_TRILI_PKG "Tpetra")
  set(PHIST_KERNEL_LIB_TPETRA 1)
elseif (PHIST_KERNEL_LIB STREQUAL "ghost")
  # we already link with ghost because of the queuing system  
  set(PHIST_KERNEL_LIB_GHOST 1)
  set(PHIST_MVECS_ROW_MAJOR 1 CACHE STRING "wether multi-vectors should be stored in row-major ordering")
  set(PHIST_SDMATS_ROW_MAJOR 0 CACHE STRING "store sdMats in row-major order (not supported by our algorithms..." )

#TODO: this is in fact matrix dependent and should be configured at run time, 
#      not compile time
  set(PHIST_USE_SELL ON CACHE BOOL "Use the GHOST SELL-C-sigma sparse matrix format")
  set(PHIST_SELL_C 1 CACHE INTEGER "chunk size C for the SELL-C-sigma format")
  set(PHIST_SELL_SIGMA 1 CACHE INTEGER "sorting scope sigma for the SELL-C-sigma format")

elseif (PHIST_KERNEL_LIB STREQUAL "builtin")
  set(PHIST_KERNEL_LIB_BUILTIN 1)
  set(PHIST_MVECS_ROW_MAJOR 1)
  option(PHIST_KERNEL_LIB_BUILTIN_PIN_THREADS "Try to pin the threads to cores using sched_setaffinity (based on simple heuristics; for Intel use env KMP_AFFINITY=\"granularity=fine,compact,1,0\"" On)
  option(PHIST_KERNEL_LIB_BUILTIN_AVX2_PREC_KERNELS "Add support for high-precision kernel operations triggered with PHIST_ROBUST_REDUCTIONS, requires AVX2,FMA (e.g. Intel Haswell)" Off)
  if(PHIST_KERNEL_LIB_BUILTIN_AVX2_PREC_KERNELS)
    message(STATUS "Activated experimental high-precision builtin kernel routines!")
    set(PHIST_HIGH_PRECISION_KERNELS 1)
    set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mavx -mavx2 -mfma")
    set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -mavx -mavx2 -mfma")
    set (CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -mavx -mavx2 -mfma")
  endif()

  if( PHIST_ENABLE_SP )
    message(FATAL_ERROR "Builtin kernel lib does not support single-precision calculations!")
  endif()
  if( NOT PHIST_HAVE_MPI )
    message(FATAL_ERROR "Builtin kernel lib requires MPI!")
  endif()
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set or not recognized" )
endif ()
set(PHIST_BENCH_LARGE_N 64000000 CACHE STRING "Dimension for benchmarking tests (e.g. large enough to make micro benchmarks memory-bound)")

# test for Belos (Trilinos library of iterative solvers).
# We recommend linking with Belos because we use its Tsqr
# interface. If Belos is not found, important functionality
# may not work, in particular the mvec_QR function is currently
# implemented via Belos by all kernel libs except "builtin"

find_package(Teuchos)

if (${Teuchos_FOUND})
  list(APPEND LOAD_TRILI_PKG "Teuchos")
  set(PHIST_HAVE_TEUCHOS 1)
  option(PHIST_USE_TEUCHOS_TIMEMONITOR "Use TimeMonitor from Trilinos package Teuchos instead of builtin one." OFF)
  list(APPEND library_list ${Teuchos_LIBRARIES})
  message(STATUS "Found Teuchos")
else()
  message(STATUS "Teuchos not found.")
endif()
option(PHIST_TIMINGS_FULL_TRACE "Try to print full call trace in the TimeMonitor/performance measurements (Warning: may generate *a lot* of output)." Off)

find_package(Kokkos)

if (${Kokkos_FOUND})
  list(APPEND LOAD_TRILI_PKG "Kokkos")
  set(PHIST_HAVE_KOKKOS 1)
  message(STATUS "Found Kokkos")
else()
  message(STATUS "Kokkos not found.")
endif()

find_package(Belos)

if (${Belos_FOUND})
  list(APPEND LOAD_TRILI_PKG "Belos")
  set(PHIST_HAVE_BELOS 1)
  message(STATUS "Found Belos")
else()
  message(STATUS "Belos not found.")
endif()

# if Anasazi is found we provide a simple interface to their
# block algorithms (e.g. Block Lanczos/Arnoldi)
find_package(Anasazi)

if (${Anasazi_FOUND})
  list(APPEND LOAD_TRILI_PKG "Anasazi")
  set(PHIST_HAVE_ANASAZI 1)
  message(STATUS "Found Anasazi")
else()
  message(STATUS "Anasazi not found.")
endif()

# even without the essex-physics libs we include certain matfuncs for the
# test cases SpinSZ and Graphene for testing purposes.
find_package(ESSEX-PHYSICS CONFIG)
if (ESSEX-PHYSICS_FOUND)
  message(STATUS "Found ESSEX-Physics ${ESSEX-PHYSICS_VERSION}")
  set(PHIST_HAVE_ESSEX_PHYSICS 1)
  include_directories(${ESSEX-PHYSICS_INCLUDE_DIRS})
else()
  message(STATUS "ESSEX-Physics not found, will use supplementary implementation")
  include_directories(${PROJECT_SOURCE_DIR}/drivers/matfuncs)
endif()

# these specific versions of type-generic drivers have to be created
if  (PHIST_KERNEL_LIB STREQUAL "epetra" OR PHIST_KERNEL_LIB STREQUAL "builtin")
set(TG_Types D)
else ()
  if (PHIST_ENABLE_SP)
    set(TG_Types S 
                 D 
                 C 
                 Z)
  else ()
    set(TG_Types D
                 Z)
  endif ()
endif ()

if (PHIST_KERNEL_LIB STREQUAL "builtin")
  set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/phist_driver_utils.c
                          ${PROJECT_SOURCE_DIR}/drivers/matfuncs/matpde.f90
                          ${PROJECT_SOURCE_DIR}/drivers/matfuncs/matpde3d.f90
                          ${PROJECT_SOURCE_DIR}/drivers/matfuncs/tritoeplitz.f90
                          ${PROJECT_SOURCE_DIR}/drivers/matfuncs/brussolator.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/phist_mpi_kernels.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/env_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/map_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/mrgrnk.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/crsmat_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/parmetis.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/colpack.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/mvec_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/nrm2_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/axpy_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/axpy_kernels_nt.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gather_scatter_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/dot_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_augmented.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sC.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_fused_sCD.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/spmvm_kernels_nt.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/spmvm_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/spmv_buff_cpy_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/kacz_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/sdmat_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/bench_kernels.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface/phist_rcp_helpers.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface/phist_BelosAdapter.cpp
                          )
  if(PHIST_KERNEL_LIB_BUILTIN_AVX2_PREC_KERNELS)
    set (PHIST_KERNEL_SOURCE ${PHIST_KERNEL_SOURCE}
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/prec_helpers.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/sdmat_kernels_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/dot_kernels_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_augmented_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sC_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_fused_sCD_prec.c)
  endif()
else()
  set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/phist_driver_utils.c
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/matpde.f90
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/matpde3d.f90
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/tritoeplitz.f90
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/brussolator.f90
                           ${PROJECT_SOURCE_DIR}/src/kernels/phist_mpi_kernels.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/builtin/bench_kernels.c
                           ${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface/phist_rcp_helpers.cpp)
endif()
if (NOT ESSEX-PHYSICS_FOUND)
  LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/drivers/matfuncs/graphene.c
                                  ${PROJECT_SOURCE_DIR}/drivers/matfuncs/SpinChainSZ.c
                                  ${PROJECT_SOURCE_DIR}/drivers/matfuncs/part_tools.c)
endif()

set (PHIST_TOOLS_SOURCE ${PROJECT_SOURCE_DIR}/src/tools/phist_tools.cpp
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_timemonitor.cpp
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_perfcheck.cpp)
set (PHIST_CORE_SOURCE ${PROJECT_SOURCE_DIR}/src/core/phist_operator.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_orthog.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_svqb.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_svrr.cpp)
set (PHIST_KRYLOV_SOURCE ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_arnoldi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedgmres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedminres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belos.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_anasazi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_carp_cg.cpp)
set (PHIST_JADA_SOURCE   ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOpts.c
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jdqr.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_schur_decomp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejada.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_transform_searchspace.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaCorrectionSolver.cpp)
set (PHIST_SUPPORT_SOURCE  ${PROJECT_SOURCE_DIR}/src/support/phist_feastCorrectionSolver.cpp
                           ${PROJECT_SOURCE_DIR}/src/support/phist_iter_op.cpp)


include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/kernels)
include_directories(${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB})
include_directories(${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface)
include_directories(${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/krylov)
include_directories(${PROJECT_SOURCE_DIR}/src/jada)
include_directories(${PROJECT_SOURCE_DIR}/src/support)
include_directories(${PROJECT_SOURCE_DIR}/drivers)

add_library(phist_tools ${PHIST_TOOLS_SOURCE})
add_library(phist_core ${PHIST_CORE_SOURCE})
add_library(phist_kernels_${PHIST_KERNEL_LIB} ${PHIST_KERNEL_SOURCE})
add_library(phist_solvers ${PHIST_KRYLOV_SOURCE} 
                          ${PHIST_JADA_SOURCE})
add_library(phist_support ${PHIST_SUPPORT_SOURCE})

LIST(APPEND library_list phist_tools)
LIST(APPEND library_list phist_core)
LIST(APPEND library_list phist_kernels_${PHIST_KERNEL_LIB})
LIST(APPEND library_list phist_solvers)

target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB}                                        phist_tools)
target_link_libraries(phist_core                             phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_solvers               phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_support phist_solvers phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)

if (ESSEX-PHYSICS_FOUND)
  target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB} ${ESSEX-PHYSICS_LIBRARIES})
endif()

if(PHIST_KERNEL_LIB STREQUAL "ghost")
  find_package(GHOST REQUIRED CONFIG)
else()
  find_package(GHOST CONFIG QUIET)
endif()
if (GHOST_FOUND)
  message(STATUS "Found GHOST ${GHOST_VERSION}")
  set(PHIST_HAVE_GHOST 1)
  include_directories(${GHOST_INCLUDE_DIRS})
  link_directories(${GHOST_LIBRARY_DIRS})
  # Moved to applications (tests and drivers)
  #target_link_libraries(phist_tools ${GHOST_LIBRARIES})
  LIST(APPEND external_libs ${GHOST_LIBRARIES})
  list(APPEND PHIST_PKG_CONFIG_REQUIRE ghost)
else()
  message(STATUS "GHOST not found, some functionality will be disabled and return -99 ('not implemented')")
endif()

message(STATUS "GHOST libs: ${GHOST_LIBRARIES}")

# load parmetis for reordering builtin
if(PHIST_KERNEL_LIB STREQUAL "builtin")
  
  find_package(ParMETIS)
  if(PARMETIS_FOUND)
    set(PHIST_HAVE_PARMETIS On)
    include_directories(${PARMETIS_INCLUDE_DIRS})
    link_directories(${PARMETIS_LIBRARY_DIRS})
    target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB} ${PARMETIS_LIBRARIES})
    message(STATUS "ParMETIS found, enabled matrix reordering!")
  else()
    message(STATUS "ParMETIS not found, no matrix reordering!")
  endif()

#  find_package(ColPack)
  if (COLPACK_DIR)
    set(COLPACK_INCLUDE_DIRS ${COLPACK_DIR}/include)
    find_library(COLPACK_LIBRARIES ColPack
                 PATHS ${COLPACK_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_DEFAULT_PATH)
    find_library(COLPACK_LIBRARY ColPack) # also check default paths, but prefer COLPACK_DIR
    #message(STATUS "COLPACK_LIBRARIES: ${COLPACK_LIBRARIES}")
    get_filename_component(COLPACK_LIBRARY_DIRS ${COLPACK_LIBRARIES} PATH)
    #message(STATUS "COLPACK_LIBRARY_DIRS: ${COLPACK_LIBRARY_DIRS}")
    #set(COLPACK_LIBRARIES ${COLPACK_DIR}/lib/libColPack.so)
    #set(COLPACK_LIBRARY_DIRS ${COLPACK_DIR}/lib)
    set(COLPACK_FOUND 1)
  endif()

  if(COLPACK_FOUND)
    set(PHIST_HAVE_COLPACK 1)
    include_directories(${COLPACK_INCLUDE_DIRS})
    link_directories(${COLPACK_LIBRARY_DIRS})
    target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB} ${COLPACK_LIBRARIES})
    message(STATUS "ColPack found, enabled matrix coloring!")
  else()
    message(STATUS "ColPack not found, no matrix coloring!")
  endif()

endif()

find_package(Likwid PATHS "$ENV{LIKWID_HOME}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" Optional)
if (${LIKWID_Found})
  set(PHIST_HAVE_LIKWID 1)
  include_directories(${LIKWID_INCLUDE_DIRS})
  target_link_libraries(phist_tools ${LIKWID_LIBRARIES})
  option(LIKWID_PERFMON "Activate function-level Likwid performance monitoring (to be used with likwid-perfctr -m)" Off)
  if (${LIKWID_PERFMON})
    message(STATUS "activated function-level Likwid performance monitoring")
  else()
    message(STATUS "disabled function-level Likwid performance monitoring")
  endif ()
else()
  message(STATUS "continuing without Likwid, you can set the environment variable LIKWID_HOME to help CMake find it")
endif()

option(PHIST_TIMEMONITOR "Gather detailed function-level timings" On)
if( PHIST_TIMEMONITOR )
  message(STATUS "activated TimeMonitor timing monitoring")
endif()
option(PHIST_TIMEMONITOR_PERLINE "Gather timing information for each line with a PHIST_CHK_* macro!" Off)
option(PHIST_PERFCHECK "Check performance of operations (which have a performance model specification)" Off)
if( PHIST_PERFCHECK )
  message(STATUS "activated PerfCheck timing monitoring")
  if( PHIST_PERFCHECK AND PHIST_TIMEMONITOR )
    message(FATAL_ERROR "You cannot use PHIST_PERFCHECK and PHIST_TIMEMONITOR at once currently!")
  endif()
  option(PHIST_PERFCHECK_REALISTIC "Switch from ideal performance predictions to more realistic ones (e.g. for views)!")
endif()


set(INCLUDE_FLAGS "")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  set(INCLUDE_FLAGS "${INCLUDE_FLAGS} -I${dir}")
endforeach()

set(PHIST_PKG_CONFIG_LIBRARIES "")
foreach(lib ${library_list})
  set(PHIST_PKG_CONFIG_LIBRARIES "${PHIST_PKG_CONFIG_LIBRARIES} -l${lib}")
endforeach()

configure_file ("${PROJECT_SOURCE_DIR}/phist.pc.in" "${PROJECT_BINARY_DIR}/phist.pc" @ONLY)
configure_file ("${PROJECT_SOURCE_DIR}/src/phist_config.h.in" "${PROJECT_BINARY_DIR}/phist_config.h")
configure_file ("${PROJECT_SOURCE_DIR}/doc/phist.doxyfile.in" "${PROJECT_BINARY_DIR}/phist.doxyfile")
include(CMakePackageConfigHelpers)
set(INCLUDE_INSTALL_DIR include)
set(LIBRARY_INSTALL_DIR lib)
configure_package_config_file("${PROJECT_SOURCE_DIR}/phistConfig.cmake.in" "${PROJECT_BINARY_DIR}/phistConfig.cmake" INSTALL_DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist
  PATH_VARS INCLUDE_INSTALL_DIR)# LIBRARY_INSTALL_DIR)
write_basic_package_version_file("${PROJECT_BINARY_DIR}/phistConfigVersion.cmake" VERSION ${PHIST_VERSION} COMPATIBILITY SameMajorVersion)

foreach(PKG ${LOAD_TRILI_PKG})
  find_package(${PKG} REQUIRED)
  if (${${PKG}_FOUND})
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    target_link_libraries(phist_tools ${${PKG}_LIBRARIES})
    target_link_libraries(phist_tools ${${PKG}_TPL_LIBRARIES})
  else()
    message(FATAL_ERROR ${PKG} " was not found")
  endif()

endforeach()

if(NOT GHOST_FOUND) #use MPI libs as found by GHOST in this case. otherwise, a clash between thread-safe and normal MPI will occur
  target_link_libraries(phist_tools ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif()
target_link_libraries(phist_tools m) #link the C-math library as it may not be automatically included in pure C programs

set(GCOVR_EXCLUDES -e src/kernels/kernels_noimpl.c
                   -e src/tools/tools.c)

# Custom properties
set(INTEGRATION_BUILD Off CACHE BOOL "Indicates an integration build")
find_program(GCOVR_SCRIPT gcovr)
find_program(CPPCHECK_EXECUTABLE cppcheck)
if (${INTEGRATION_BUILD})
    set(GCOVR_OPTIONS -r src ${GCOVR_EXCLUDES} --xml -o coverage.xml CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS --force --enable=all --suppress=unusedFunction --xml --xml-version=2 -q CACHE_STRING "cppcheck options")
    set(CPPCHECK_OPTIONS_POST 2> cppcheck.xml)
    set(CTEST_OPTIONS --force-new-ctest-process --verbose CACHE STRING "ctest options")
else()
    set(GCOVR_OPTIONS -r src ${GCOVR_EXCLUDES} CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS --force --enable=all --suppress=unusedFunction CACHE STRING "cppcheck options")
    set(CPPCHECK_OPTIONS_POST)
    set(CTEST_OPTIONS --force-new-ctest-process --output-on-failure CACHE STRING "ctest options")
endif()

# Define includes
include_directories(lib/googletest-1.6 src/ src/tools/ src/kernels/)

# Define sources
set(SOURCES 
  ${PHIST_KERNEL_SOURCE}
  ${PHIST_TOOLS_SOURCE}
  ${PHIST_CORE_SOURCE}
  ${PHIST_KRYLOV_SOURCE}
  ${PHIST_JADA_SOURCE}
  ${PHIST_SUPPORT_SOURCE}
)

set(TEST_SOURCES
    lib/googletest-1.6/gmock-gtest-all.cc
    test/mainTest.cpp
    test/tools/MatrixIO.cpp)
    
file(GLOB TestMatricesMM ${CMAKE_SOURCE_DIR}/test/matrices/*.mm)
file(GLOB TestMatricesBIN ${CMAKE_SOURCE_DIR}/test/matrices/*.bin)

foreach (mat ${TestMatricesMM})
  configure_file(${mat} ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
  list(APPEND TEST_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${mat_name})
endforeach()
foreach (mat ${TestMatricesBIN})
  configure_file(${mat} ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
  list(APPEND TEST_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${mat_name})
endforeach()

#######################
# Define executables  #
#######################

# Directly link all external libraries into tools as all other
# phist libraries may depend on it and the *user* cannot now
# what libraries to link. (melven)
# Additionally, this is required to make the python interface work...
message(STATUS "External libraries linked: ${external_libs}")
target_link_libraries(phist_tools ${external_libs})

# type generic driver routines are compiled from .c or .cpp files
# which we get from .in files by preprocessing
file(GLOB TG_C_Drivers ${CMAKE_SOURCE_DIR}/drivers/*.c.in )
file(GLOB TG_CXX_Drivers ${CMAKE_SOURCE_DIR}/drivers/*.cpp.in )

foreach (infile ${TG_C_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.c)
    add_executable(phist_${tc}${filebase} ${tc}${filebase}.c)
    target_link_libraries(phist_${tc}${filebase} ${library_list})
    install(TARGETS phist_${tc}${filebase} RUNTIME DESTINATION bin OPTIONAL)
    # added
    LIST(APPEND driver_list phist_${tc}${filebase})
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(phist_${tc}${filebase} ${tc}${filebase}.cpp)
    target_link_libraries(phist_${tc}${filebase} ${library_list})
    install(TARGETS phist_${tc}${filebase} RUNTIME DESTINATION bin OPTIONAL)
    LIST(APPEND driver_list phist_${tc}${filebase})
  endforeach()
endforeach()

# example drivers using matfuncs

add_executable(phist_carp_cg ${PROJECT_SOURCE_DIR}/drivers/carp_cg.c)
target_link_libraries(phist_carp_cg ${library_list})

# these drivers need the feastCorrectionSolver and/or the iter_op
target_link_libraries(phist_carp_cg phist_support)
foreach (tc ${TG_Types})
  target_link_libraries(phist_${tc}krylov_schur phist_support)
endforeach()


install(TARGETS phist_carp_cg RUNTIME DESTINATION bin OPTIONAL)

LIST(APPEND driver_list phist_carp_cg)

if (${BUILD_SHARED_LIBS})
  install(TARGETS phist_kernels_${PHIST_KERNEL_LIB} EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_tools EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_core EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_solvers EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
else()
  install(TARGETS phist_kernels_${PHIST_KERNEL_LIB} EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_tools EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_core EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_solvers EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
endif()
install(FILES ${PROJECT_BINARY_DIR}/phistConfig.cmake        DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist)
#export(EXPORT phistLibraries FILE "${PROJECT_BINARY_DIR}/phistLibraries.cmake")
install(EXPORT phistLibraries FILE phistLibraries.cmake DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist)
install(FILES ${PROJECT_BINARY_DIR}/phistConfigVersion.cmake DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist)
install(FILES ${PROJECT_BINARY_DIR}/phist.pc DESTINATION ${LIBRARY_INSTALL_DIR}/pkgconfig)
install(FILES ${PROJECT_BINARY_DIR}/phist_config.h DESTINATION ${INCLUDE_INSTALL_DIR})

file(GLOB PHIST_HEADERS ${CMAKE_SOURCE_DIR}/src/*/*.h 
                        ${CMAKE_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/phist_typedefs.h 
                        ${CMAKE_SOURCE_DIR}/src/tools/*.hpp
                        )
install(FILES ${PHIST_HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR})

file(GLOB PYTHON_HEADERS ${CMAKE_SOURCE_DIR}/python/*.py)
install(FILES ${PYTHON_HEADERS} DESTINATION ${LIBRARY_INSTALL_DIR}/python)

# test drivers
if( "${CMAKE_BUILD_TYPE}" MATCHES "Release" )
  set(PHIST_TESTS task
                  kernels
                  core
                  krylov
                  jada
                  bench)
else()
  set(PHIST_TESTS task
                  kernels
                  core
                  krylov
                  jada)
endif()

foreach(label ${PHIST_TESTS})
  string(TOUPPER ${label} label_uc)

  file(GLOB ${label_uc}_TEST_SOURCES ${CMAKE_SOURCE_DIR}/test/${label}/*.cpp)
  list(APPEND ${label_uc}_TEST_SOURCES ${TEST_SOURCES})

  add_executable(${PRODUCT_NAME}-${label}-test ${${label_uc}_TEST_SOURCES} ${TEST_DEPENDS})
  set_target_properties(${PRODUCT_NAME}-${label}-test PROPERTIES COMPILE_DEFINITIONS "PHIST_TESTSUITE=\"${label}\"")
  # compile tests without optimization
  set_target_properties(${PRODUCT_NAME}-${label}-test PROPERTIES COMPILE_FLAGS "-O0")
  target_link_libraries(${PRODUCT_NAME}-${label}-test ${library_list})
  target_link_libraries(${PRODUCT_NAME}-${label}-test pthread)
  LIST(APPEND test_list ${PRODUCT_NAME}-${label}-test)

endforeach()

# some example drivers which haven't made their way into the drivers/ directory
file(GLOB C_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.c )
file(GLOB CXX_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.cpp )

foreach (infile ${C_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  add_executable(phist_${filebase} ${infile})
  target_link_libraries(phist_${filebase} ${library_list})
  LIST(APPEND test_list phist_${filebase})
endforeach()

foreach (infile ${CXX_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  add_executable(phist_${filebase} ${infile})
  target_link_libraries(phist_${filebase} ${library_list})
  LIST(APPEND test_list phist_${filebase})
endforeach()

# type generic examples routines are compiled from .c or .cpp files
# which we get from .in files by preprocessing
file(GLOB TG_C_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.c.in )
file(GLOB TG_CXX_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.cpp.in )

foreach (infile ${TG_C_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.c)
    add_executable(${tc}${filebase} ${tc}${filebase}.c)
    target_link_libraries(${tc}${filebase} ${library_list})
    LIST(APPEND test_list ${tc}${filebase})
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(${tc}${filebase} ${tc}${filebase}.cpp)
    target_link_libraries(${tc}${filebase} ${library_list})
    LIST(APPEND test_list ${tc}${filebase})
  endforeach()
endforeach()

#

# Activate testing
enable_testing()

if (${PHIST_ENABLE_MPI})
  set(PHIST_MPI_TEST_NPROCS 1 2 4 5 6)
else()
  set(PHIST_MPI_TEST_NPROCS "")
endif()

if (PHIST_ENABLE_COMPLEX_TESTS)
  set(GTEST_FILTER "")
else()
  set(GTEST_FILTER "--gtest_filter=-C*:Z*")
endif()

foreach(np ${PHIST_MPI_TEST_NPROCS})
  math(EXPR TEST_OMP_NUM_THREADS "6/${np}")
  foreach(label ${PHIST_TESTS})
    set(TEST_COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} env LSAN_OPTIONS=\"exitcode=0\" env PHIST_NUM_THREADS=${TEST_OMP_NUM_THREADS} env OMP_NUM_THREADS=${TEST_OMP_NUM_THREADS} env OMP_SCHEDULE=guided,256
    ./${PRODUCT_NAME}-${label}-test --gtest_death_test_style=threadsafe --gtest_output=xml:test_output_np${np}/${label}/xunit.xml ${GTEST_FILTER})
    add_test(NAME test-${label}_np${np} COMMAND ${TEST_COMMAND})
    set_tests_properties(test-${label}_np${np} PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*PASSED[. \n]*")
    set_tests_properties(test-${label}_np${np} PROPERTIES FAIL_REGULAR_EXPRESSION "[. \n]*FAILED[. \n]*")
    set_tests_properties(test-${label}_np${np} PROPERTIES TIMEOUT 500)
  endforeach()
  # the task tests should be very fast (no heavy computation) - but may deadlock in case of errors!
  set_tests_properties(test-task_np${np} PROPERTIES TIMEOUT 10)
endforeach()

# Test python interface if python is available
find_package(PythonInterp)
if(PYTHONINTERP_FOUND)
  foreach(pyh ${PYTHON_HEADERS})
    get_filename_component(pyh_name ${pyh} NAME)
    add_test(NAME test-${pyh_name} COMMAND ${PYTHON_EXECUTABLE} ${pyh})
    set_tests_properties(test-${pyh_name} PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH}")
    if( GCC_SANITIZE )
      # we need to LD_PRELOAD the respective lib<.>san.so library
      string(SUBSTRING ${GCC_SANITIZE} 0 1 GCC_SANITIZE_FIRST_CHAR)
      set(GCC_SANITIZE_LIB "lib${GCC_SANITIZE_FIRST_CHAR}san.so")
      # work around cmake bug (specifying multiple set_tests_properties(...ENVIRONMENT...)
      set_property(TEST test-${pyh_name} APPEND PROPERTY ENVIRONMENT "LD_PRELOAD=${GCC_SANITIZE_LIB}:$ENV{LD_PRELOAD}")
      set_property(TEST test-${pyh_name} APPEND PROPERTY ENVIRONMENT "LSAN_OPTIONS=\"exitcode=0\"")
    endif()
  endforeach()
endif()

# Adding concurrency tests with valgrind-helgrind.
# This was done by Timo for his bachelor thesis, currently
# the whole tasking thing is not working and the use/portability
# of automated valgrind tests is questionable, so we skip this test for now.
#add_test(concurrency-test: ${CMAKE_SOURCE_DIR}/test/multiple_helgrind_run.sh ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR})

# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS})

# Activate gcov code generation
if (CMAKE_COMPILER_IS_GNUCXX AND "${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    message (STATUS "Coverage support enabled for targets: ${ARG_TARGETS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
endif()

# Custom targets: doc, cppcheck, coverage, audit
# Adding optinal target "doc" to generate Doxygen documentation
find_package(Doxygen 1.7.0)
if (DOXYGEN_FOUND)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/phist.doxyfile WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# Add optional target "coverage" to genereate code coverage report (console / XML)
if (GCOVR_SCRIPT)
    add_custom_target(coverage ${GCOVR_SCRIPT} ${GCOVR_OPTIONS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# Adding optional "cppcheck" target
if (CPPCHECK_EXECUTABLE)
    # get all C and C++ sources
    set(C_SOURCES)
    foreach(source ${SOURCES})
      string(REGEX MATCH ".*[.]cp?p?" c_source ${source})
      if(c_source)
        list(APPEND C_SOURCES ${c_source})
      endif()
    endforeach()
    add_custom_target(cppcheck ${CPPCHECK_EXECUTABLE} ${CPPCHECK_OPTIONS} ${INCLUDE_FLAGS} ${C_SOURCES} ${CPPCHECK_OPTIONS_POST})
endif()

# Audit target: Combines check, cppcheck, coverage
add_custom_target(audit DEPENDS coverage cppcheck)

# Define the distribution package
set(CPACK_PACKAGE_VERSION "${PHIST_VERSION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_C_COMPILER_ID}-${CMAKE_C_COMPILER_VERSION}")
set(CPACK_GENERATOR TGZ)

install(FILES "${CMAKE_SOURCE_DIR}/LICENSE.txt" DESTINATION share/doc/phist)
install(FILES "${CMAKE_SOURCE_DIR}/README.md"   DESTINATION share/doc/phist)
if (DOXYGEN_FOUND)
  install(DIRECTORY "${CMAKE_SOURCE_DIR}/doc/html" DESTINATION share/doc/phist OPTIONAL)
endif()
include(CPack)

# Custom targets
add_custom_target(libs DEPENDS ${library_list})
add_custom_target(drivers DEPENDS ${driver_list})
add_custom_target(tests DEPENDS ${test_list})
