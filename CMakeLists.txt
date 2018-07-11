
cmake_minimum_required(VERSION 3.8.0)

# some guidelines for extending/modifying the PHIST CMakeLists.txt file. 
# Things may not be handled this way everywhere yet, but you should keep these
# points in mind.
#
# - anything prefixed by PHIST_ENABLE or PHIST_USE will enter the CMake Cache but not the config file,
#   so it won't be directly visible in the source code
# - anything prefixed by PHIST_HAVE can be seen in the phist_config.h file
# - other things should *not* go into the cache unless they come from another library via find_package.
#
# - try to mark things that are very detailed (like details on the timing behavior etc.) as ADVANCED so that
#   the CMake cache seen when first using ccmake stays reasonable

#######################################################################
# PROJECT SETTINGS
#######################################################################

# Define some custom variables
project(phist C CXX Fortran)

# here an odd version minor (e.g. 1.1.x, 1.3.x etc.) indicates a development version
# with unstable functionality, an even version minor (1.0.x, 1.2.x etc.) a stable version.
# Stable versions will return *iflag=-99 (not implemented) if a function is called that is
# still under development, whereas the next development version may or may not do what
# the function claims to do, depending on the status of the repository.
set(PHIST_VERSION_MAJOR 1)
set(PHIST_VERSION_MINOR 7)
set(PHIST_VERSION_PATCH 1)
set(PHIST_VERSION ${PHIST_VERSION_MAJOR}.${PHIST_VERSION_MINOR}.${PHIST_VERSION_PATCH})
math(EXPR PHIST_VERSION_INT "( ( ${PHIST_VERSION_MAJOR} * 100 ) + ${PHIST_VERSION_MINOR} ) * 100 + ${PHIST_VERSION_PATCH}")

set(PRODUCT_NAME phist-${PHIST_VERSION})

execute_process(COMMAND git rev-parse --short --sq HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE rev)
if (rev)
    STRING(REPLACE "'" "\"" PHIST_GIT_REVISION ${rev})
else()
    SET(PHIST_GIT_REVISION "\"${PHIST_VERSION}\"")
endif()

##########################################################################################
# GENERAL CMAKE SETTINGS
##########################################################################################

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# install target should not depend on 'make all': you can do e.g. 'make libs; make install'
SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Release, RelWithDebug, Debug" FORCE)

set(BUILD_SHARED_LIBS ON CACHE BOOL "Set to ON/OFF for shared/static library build.")

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##########################################################################################
# OPTIONS FOR COMPILE AND RUN TIME
##########################################################################################

# setup which kernel library is being used.
set(PHIST_KERNEL_LIB "builtin" CACHE STRING "Select underlying kernel libraries for phist (builtin|ghost|epetra|tpetra|magma|petsc|eigen)")
message(STATUS "PHIST kernel library: " ${PHIST_KERNEL_LIB})

include(cmake/PhistOptions.cmake)

########################################################################
# XSDK OPTIONS
########################################################################

include(cmake/XsdkOptions.cmake)

########################################################################
# try to enable CCache as a 'compiler launcher'
########################################################################

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND AND PHIST_USE_CCACHE)
    message(STATUS "enabling CCache to speed-up repeated builds.")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif()
endif()

########################################################################
# Lists of libraries, executables, tests
########################################################################
set(library_list)
set(driver_list)
set(test_list)

########################################################################
# List of external libraries for executables and tests to link, e.g.
# GHOST, LAPACK.
########################################################################
set(external_libs)
set(PHIST_PKG_CONFIG_REQUIRE)
set(PHIST_PKG_CONFIG_CFLAGS)

########################################################################
# Deal with BLAS, LAPACK and LAPACKe
#
# As required by the xSDK, the user can explicitly set the following variables and we will respect them:
#
# TPL_LAPACKE_LIBRARIES (which should include the lapacke, lapack and blas libs with full names (no -L and -l options)
# TPL_LAPACKE_INCLUDE_DIR
#
# and to conform to standard flags, we also allow setting
#
# TPL_LAPACK_LIBRARIES (which should include BLAS libraries)
#
# and will attempt to find lapacke ourselves.
#
# If none of the above are given, we try to find the libraries ourselves.
# For the Intel compiler we set the linker flag "-mkl=sequential".
# Otherwise, we trust the CMake FindLAPACK module.
# In order to use a specific blas/lapack (but you don't know the exact paths), 
# you can set the variable BLA_VENDOR used in that module, e.g.:
# BLA_VENDOR="Intel10_64lp" (note that you have to be rather specific here
# due to the FindBLAS module!).
# There is no FindLAPACKe module as of CMake 3.8, so I added one in phist/cmake/
# If MKL is used, skip lapacke (because it is included in MKL) but explicitly add 
# the include directory.
########################################################################

set(TPL_LAPACKE_LIBRARIES "${TPL_LAPACKE_LIBRARIES}" CACHE STRING "full library names of lapacke, lapack and blas" FORCE)
set(TPL_LAPACK_LIBRARIES "${TPL_LAPACK_LIBRARIES}" CACHE STRING "full library names of lapack and blas" FORCE)
set(TPL_LAPACKE_INCLUDE_DIRS "${TPL_LAPACKE_INCLUDE_DIRS}" CACHE STRING "path where to find lapacke.h or mkl_lapacke.h header" FORCE)

if (TPL_LAPACKE_LIBRARIES)
  list(GET TPL_LAPACKE_LIBRARIES 0 LAPACKE_LIBRARY)
  get_filename_component(LAPACK_LIBDIR ${LAPACKE_LIBRARY} DIRECTORY)
elseif(TPL_LAPACK_LIBRARIES)
  list(GET TPL_LAPACK_LIBRARIES 0 LAPACK_LIBRARY)
  get_filename_component(LAPACK_LIBDIR ${LAPACK_LIBRARY} DIRECTORY)
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(TPL_LAPACKE_LIBRARIES "MKL_VIA_FLAG")
  set(TPL_LAPACKE_INCLUDE_DIRS "MKL_VIA_FLAG")
  set(MKL_FLAG "-mkl=sequential")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKL_FLAG}")
else()
  # try to find lapack and blas using the standard CMake module.
  # It searches for various vendors including Intel MKL, GOTO BLAS etc.
  # To specify the lapack library to search for, the user can set the
  # BLA_VENDOR variable either via -D or as environment variable.
  find_package(LAPACK REQUIRED QUIET)
  link_libraries(${LAPACK_LIBRARIES})
  list(GET LAPACK_LIBRARIES 0 LAPACK_LIBRARY)
  get_filename_component(LAPACK_LIBDIR ${LAPACK_LIBRARY} DIRECTORY)
  set(TPL_LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE STRING "lapack library and dependencies" FORCE)
endif()

if (NOT TPL_LAPACKE_LIBRARIES)
  if (NOT TPL_LAPACK_LIBRARIES)
    message(FATAL_ERROR "did not find lapack via TPL_LAPACK_LIBRARIES, TPL_LAPACK_LIBRARIES or find_package. Aborting.")
  elseif (TPL_LAPACK_LIBRARIES MATCHES "mkl")
    # note: lapacke is included in the mkl library
    set(TPL_LAPACKE_LIBRARIES ${TPL_LAPACK_LIBRARIES} CACHE STRING "full library names of lapack and blas." FORCE)
  else()
    find_library (LAPACKE_LIBRARY lapacke ${LAPACKE_LIBDIR})
    if (LAPACKE_LIBRARY MATCHES "LAPACKE_LIBRARY-NOTFOUND")
      message(FATAL_ERROR "Did not find the LAPACKe library based on the LAPACK_LIBRARIES you provided, please set "
                          "LAPACKE_LIBRARIES and LAPACKE_INCLUDE_DIRS instead (including lapacke,lapack and blas library "
                          "names/inc dirs) or leave it to CMake to locate BLAS, LAPACK and LAPACKe.")
    else()
      message(STATUS "found LAPACKe library: '${LAPACKE_LIBRARY}'")
      set(TPL_LAPACKE_LIBRARIES ${LAPACKE_LIBRARY} ${TPL_LAPACK_LIBRARIES} CACHE STRING "LAPACKe library and dependencies" FORCE)
    endif()
  endif()
endif()
message(STATUS "TPL_LAPACKE_LIBRARIES has final value '${TPL_LAPACKE_LIBRARIES}'")
if (NOT TPL_LAPACKE_LIBRARIES)
    message(FATAL_ERROR "did not find LAPACKe via TPL_LAPACK_LIBRARIES or manually. Aborting.")
endif()

if (TPL_LAPACKE_LIBRARIES MATCHES "mkl" OR TPL_LAPACKE_LIBRARIES MATCHES "MKL_VIA_FLAG")
  message(STATUS "found Intel(R) MKL")
  set(PHIST_HAVE_MKL 1)
endif()
if (NOT TPL_LAPACKE_LIBRARIES MATCHES "MKL_VIA_FLAG")
  link_libraries(${TPL_LAPACKE_LIBRARIES})
endif()

# find the lapacke.h header (for MKL it's called mkl_lapacke.h). We use the path to the lapack lib as a hint.
if (NOT TPL_LAPACKE_INCLUDE_DIRS)
  find_path (LAPACKE_INCLUDE_DIR NAMES lapacke.h mkl_lapacke.h 
                                 HINTS ${LAPACK_LIBDIR}
                                 PATH_SUFFIXES "lapacke")
  if (LAPACKE_INCLUDE_DIR MATCHES "LAPACKE_INCLUDE_DIR-NOTFOUND")
    message(WARNING "did not find lapacke header, if the build fails you may need to set TPL_LAPACKE_INCLUDE_DIRS explicitly.")
  else()
    message(STATUS "found LAPACKe header in '${LAPACKE_INCLUDE_DIR}'")
    set(TPL_LAPACKE_INCLUDE_DIRS ${LAPACKE_INCLUDE_DIR})
  endif()
endif()
if (NOT LAPACKE_INCLUDE_DIRS MATCHES "MKL_VIA_FLAG")
  include_directories(${TPL_LAPACKE_INCLUDE_DIRS})
endif()

########################################################################
# We need Python e.g. for generating the Fortran and C++ bindings      #
########################################################################
find_package(PythonInterp 3)
if (NOT PythonInterp_FOUND)
  message(WARNING "Did not find a Python3 interpreter, will disable Fortran, C++ and Python bindings and skip some tests.")
  set(XSDK_ENABLE_Fortran OFF)
  set(XSDK_ENABLE_CXX OFF)
  set(DISABLE_PYTHON_TESTS 1)
endif()

##############################################################################
# compiler-dependent flags 
##############################################################################

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")

  message(STATUS "using GNU compiler flags")
  include(cmake/SetCompilerFlagsGNU.cmake)
  
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")

  message(STATUS "using INTEL compiler flags")
  include(cmake/SetCompilerFlagsINTEL.cmake)

elseif (CMAKE_C_COMPILER_ID STREQUAL "PGI")

  message(STATUS "using PGI compiler flags")
  include(cmake/SetCompilerFlagsPGI.cmake)

elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")

  message(STATUS "using Clang (LLVM) compiler flags")
  include(cmake/SetCompilerFlagsCLANG.cmake)

else()

  message ("compiler ID: " ${CMAKE_C_COMPILER_ID})
  message ("No custom compiler flags are known.")

endif ()

# detect if certain C++ and Fortran features are supported
include(cmake/CheckLanguageFeatures.cmake)

if(PHIST_HAVE_CXX11_LAMBDAS OR PHIST_HAVE_CXX11_THREADLOCAL)
  set(PHIST_PKG_CONFIG_CFLAGS "${PHIST_PKG_CONFIG_CFLAGS} -std=c++11")
endif()

# check for SSE, AVX, AVX2 and AVX512
include(cmake/CheckCPUFeatures.cmake)

if (${PHIST_ENABLE_MPI})
  set(PHIST_HAVE_MPI 1)
#TODO - what is the best way to do this portably?
  find_package(MPI QUIET)
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
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  endif()
  set(MPIEXEC "mpirun" CACHE STRING "")
  set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")
  set(MPIEXEC_ENVIRONMENT "OMPI_MCA_hwloc_base_binding_policy=none" CACHE STRING "")
else()
  set(PHIST_HAVE_MPI 0)
endif()

set(REQUIRE_TRILI_PKG,"")
set(OPTIONAL_TRILI_PKG,"")

if (PHIST_USE_SOLVER_TPLS)
  list(APPEND OPTIONAL_TRILI_PKG "Teuchos")
  list(APPEND OPTIONAL_TRILI_PKG "Anasazi")
  list(APPEND OPTIONAL_TRILI_PKG "Belos")
endif()

if  (PHIST_KERNEL_LIB STREQUAL "epetra")
  list(APPEND REQUIRE_TRILI_PKG "Epetra")
  list(APPEND REQUIRE_TRILI_PKG "EpetraExt")
  list(APPEND REQUIRE_TRILI_PKG "Teuchos")
  set(PHIST_KERNEL_LIB_EPETRA 1)

  list(APPEND OPTIONAL_TRILI_PKG "Kokkos")
  if (PHIST_USE_GRAPH_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Isorropia")
    list(APPEND OPTIONAL_TRILI_PKG "Zoltan")
  endif()
  if (PHIST_USE_PRECON_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Ifpack")
    list(APPEND OPTIONAL_TRILI_PKG "ML")
    list(APPEND OPTIONAL_TRILI_PKG "MueLU")
  endif()
elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")
  list(APPEND REQUIRE_TRILI_PKG "Tpetra")
  list(APPEND REQUIRE_TRILI_PKG "Kokkos")
  list(APPEND REQUIRE_TRILI_PKG "Teuchos")
  set(PHIST_KERNEL_LIB_TPETRA 1)
  if (PHIST_USE_PRECON_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Ifpack2")
    list(APPEND OPTIONAL_TRILI_PKG "Amesos2")
    list(APPEND OPTIONAL_TRILI_PKG "MueLU")
  endif()
elseif (PHIST_KERNEL_LIB STREQUAL "ghost")
  # we already link with ghost because of the queuing system  
  set(PHIST_KERNEL_LIB_GHOST 1)
  set(PHIST_USE_GHOST ON)
  set(PHIST_MVECS_ROW_MAJOR 1 CACHE STRING "wether multi-vectors should be stored in row-major ordering")
  set(PHIST_USE_GHOST_TASKS 0 CACHE STRING "execute all kernels as separate GHOST tasks, required for using PHISTs tasking features.")

  set(PHIST_SELL_C -1 CACHE STRING "chunk size C for the SELL-C-sigma format (-1: automatic choice)")
  set(PHIST_SELL_SIGMA -1 CACHE STRING "sorting scope sigma for the SELL-C-sigma format (-1: automatic choice)")

  set(PHIST_DEFAULT_SPMV_MODE 0 CACHE STRING "default GHOST spmv mode (0: VECTOR, 8: OVERLAP, 16: TASK)")

  if (PHIST_USE_GRAPH_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Zoltan")
  endif()
  list(APPEND OPTIONAL_TRILI_PKG "Kokkos")
  list(APPEND OPTIONAL_TRILI_PKG "Teuchos")

elseif (PHIST_KERNEL_LIB STREQUAL "builtin")
  set(PHIST_KERNEL_LIB_BUILTIN 1)
  set(PHIST_MVECS_ROW_MAJOR 1)
  if(PHIST_HAVE_AVX2)
    option(PHIST_HIGH_PRECISION_KERNELS "Activate experimental high-precision builtin kernel routines" Off)
    option(PHIST_HIGH_PRECISION_KERNELS_FORCE "Always force the usage of high precision reductions even if this requires copying data and is slow!" Off)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mavx2 -mfma")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mfma")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx -mavx2 -mfma")
  endif()

  if( NOT PHIST_HAVE_MPI )
    message(FATAL_ERROR "Builtin kernel lib requires MPI!")
  endif()
elseif (PHIST_KERNEL_LIB STREQUAL "magma")
  set(PHIST_KERNEL_LIB_MAGMA 1)
  set(PHIST_MVECS_ROW_MAJOR 0)
elseif (PHIST_KERNEL_LIB STREQUAL "petsc")
  set(PHIST_KERNEL_LIB_PETSC 1)
  set(PHIST_MVECS_ROW_MAJOR 0)
elseif (PHIST_KERNEL_LIB STREQUAL "eigen")
  set(PHIST_KERNEL_LIB_EIGEN 1)
  set(PHIST_MVECS_ROW_MAJOR 0 CACHE STRING "wether multi-vectors should be stored in row-major ordering")
  if( PHIST_HAVE_MPI )
    message(FATAL_ERROR "The eigen kernel lib has no MPI support!")
  endif()
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set or not recognized" )
endif ()
set(PHIST_BENCH_LARGE_N 134217728 CACHE STRING "Dimension for benchmarking tests (e.g. large enough to make micro benchmarks memory-bound). Should be a multiple of 8.")

option(PHIST_TIMINGS_FULL_TRACE "Try to print full call trace in the TimeMonitor/performance measurements (Warning: may generate *a lot* of output)." Off)

if (NOT DEFINED PHIST_FORCE_32BIT_GIDX_DEFAULT)

#############################################################################
# use 32-bit integers as global indices with Epetra by default. This is     #
# because e.g. Isorropia/Zoltan and our own HYMLS require it. For larger    #
# runs with Epetra it has to be switched off before compiling PHIST.        #
#############################################################################
  if( "${PHIST_KERNEL_LIB}" MATCHES "epetra|petsc|eigen" )
    set(PHIST_FORCE_32BIT_GIDX_DEFAULT ON)
  else()
    set(PHIST_FORCE_32BIT_GIDX_DEFAULT OFF)
  endif()
endif()

option(PHIST_FORCE_32BIT_GIDX "Use 32-bit global indices for compatiblity with other libs" 
${PHIST_FORCE_32BIT_GIDX_DEFAULT})

##############################################################################
# enable single precision and complex valued calculation if supported
##############################################################################
if( "${PHIST_KERNEL_LIB}" MATCHES "tpetra|ghost|magma|petsc|eigen|tpetra12" )
  set(PHIST_SUPPORT_SP ON)
  set(PHIST_SUPPORT_COMPLEX ON)
else()
  set(PHIST_SUPPORT_SP OFF)
  set(PHIST_SUPPORT_COMPLEX OFF)
endif()

option(PHIST_ENABLE_SP "Enable single precision functions and drivers" ${PHIST_SUPPORT_SP})
option(PHIST_ENABLE_COMPLEX "Enable complex valued functions and drivers" ${PHIST_SUPPORT_COMPLEX})

if (PHIST_ENABLE_SP)
  if(NOT PHIST_SUPPORT_SP)
    message(WARNING "Enabled PHIST_ENABLE_SP, but ${PHIST_KERNEL_LIB} kernel lib does not support single-precision calculations!")
  endif()
  set(PHIST_HAVE_SP 1)
endif ()

if (PHIST_ENABLE_COMPLEX)
  if(NOT PHIST_SUPPORT_COMPLEX)
    message(WARNING "Enabled PHIST_ENABLE_COMPLEX, but ${PHIST_KERNEL_LIB} kernel lib does not support complex-valued calculations!")
  endif()
  set(PHIST_HAVE_CMPLX 1)
endif()

########################################################
# configure scamac subtree (scalable matrix generator) #
########################################################

# some matfuncs (for row-wise scalable generation of test matrices) are
# included in phist. In particular we now have the 'Scalable Matrix Collection' (scamac)
# as a subtree. The stand-alone repo is at https://bitbucket.org/essex/MatrixCollection.
include_directories(${PROJECT_SOURCE_DIR}/drivers/matfuncs)
include_directories(${PROJECT_SOURCE_DIR}/drivers/matfuncs/scamac)
include_directories(${CMAKE_BINARY_DIR}/drivers/matfuncs/scamac)
set(SCAMAC_BUILD_MWE OFF CACHE BOOL "")
# to avoid conflicts with e.g. ghost (which also uses cblas), disable blas and lapack in the scamac within phist.
# they are used for allowing to compute some eigenvalues etc., but we are only interested in the matrix generator
# itself.
set(SCAMAC_USE_BLAS OFF CACHE BOOL "")
set(SCAMAC_USE_LAPACK OFF CACHE BOOL "")
if (PHIST_FORCE_32BIT_GIDX)
  set(SCAMAC_USE_64 OFF CACHE BOOL "")
else()
  set(SCAMAC_USE_64 ON CACHE BOOL "")
endif()
add_subdirectory(drivers/matfuncs/scamac)


########################################################
# define source files and libraries                    #
########################################################

# these specific versions of type-generic drivers have to be created
if (PHIST_ENABLE_SP)
  if (PHIST_ENABLE_COMPLEX)
    set(TG_Types S D C Z)
  else ()
    set(TG_Types S D)
  endif ()
else ()
  if (PHIST_ENABLE_COMPLEX)
    set(TG_Types D Z)
  else ()
    set(TG_Types D)
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
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_add_sD.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/spmvm_kernels_nt.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/spmvm_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/spmv_buff_cpy_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/kacz_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/sdmat_module.f90
                          )
  if(PHIST_HIGH_PRECISION_KERNELS)
    set (PHIST_KERNEL_SOURCE ${PHIST_KERNEL_SOURCE}
                          ${PROJECT_SOURCE_DIR}/src/tools/phist_prec_helpers.c
                          ${PROJECT_SOURCE_DIR}/src/core/DsdFact_prec_kernels.c
                          ${PROJECT_SOURCE_DIR}/src/core/mlapack_wrapper.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/sdmat_kernels_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/dot_kernels_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_augmented_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sB_add_sD_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_sC_prec.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/gemm_kernels_fused_sCD_prec.c)
  endif()
else()
  set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernelsS.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernelsC.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernelsD.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernelsZ.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/phist_driver_utils.c
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/matpde.f90
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/matpde3d.f90
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/tritoeplitz.f90
                           ${PROJECT_SOURCE_DIR}/drivers/matfuncs/brussolator.f90
                           ${PROJECT_SOURCE_DIR}/src/kernels/phist_mpi_kernels.cpp)
endif()

# these matfuncs have always been in phist and remain there for backwards compatibility,
# even though the scamac can also generate them.
LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/drivers/matfuncs/graphene.c
                                  ${PROJECT_SOURCE_DIR}/drivers/matfuncs/SpinChainSZ.c
                                  ${PROJECT_SOURCE_DIR}/drivers/matfuncs/part_tools.c)

# repartitioning with Isorropia/Zoltan
if (PHIST_KERNEL_LIB STREQUAL "epetra")
  LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/epetra/hg_repart.cpp)
endif()

# some `derived' kernel functions, e.g. mvec_my_length is implemented via mvec_get_map/map_get_local_lenth etc.
# and common things like benchmarks
LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/common/kernels_common_impl.cpp
                                ${PROJECT_SOURCE_DIR}/src/kernels/common/stream_bench.c)

if (PHIST_USE_CRAFT)
  set(CRAFT_DIR "${PROJECT_SOURCE_DIR}/src/craft/") 
  set(CRAFT_BUILD_GHOST OFF)
  set(CRAFT_BUILD_PHIST ON)
  include(${CRAFT_DIR}/CRAFT_VARS.txt) 
  add_definitions(-DPHIST_CP)
  add_library(craft SHARED ${CRAFT_SOURCE}) 
endif()

set (PHIST_TOOLS_SOURCE ${PROJECT_SOURCE_DIR}/src/tools/phist_tools.cpp
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_random.c
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_timemonitor.cpp
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_perfcheck.cpp
                        )
set (PHIST_CORE_SOURCE ${PROJECT_SOURCE_DIR}/src/core/phist_operator.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_normF.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_sdFact.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_orthog.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_svqb.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_chol_QR.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_orthogrr.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_orthogrrfused.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_svrr.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_schur_decomp.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_transform_searchspace.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_kpm.cpp
                       )

if (PHIST_KERNEL_LIB STREQUAL "epetra")
  list (APPEND PHIST_CORE_SOURCE ${PROJECT_SOURCE_DIR}/src/core/phist_operator_epetra.cpp)
endif()

set (PHIST_KRYLOV_SOURCE ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_arnoldi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_lanczos.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedgmres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedqmr.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedbicgstab.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedminres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belosS.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belosC.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belosD.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belosZ.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_anasaziS.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_anasaziC.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_anasaziD.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_anasaziZ.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_carp_cg.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_iter_op.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedpcg.cpp
                         )
if(PHIST_USE_CRAFT)
        set (PHIST_KRYLOV_SOURCE ${PHIST_KRYLOV_SOURCE}
                                                 ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_lanczos_ft.cpp
                                                                                                 )
endif()
set (PHIST_JADA_SOURCE   ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOpts.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejadaS.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejadaC.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejadaD.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejadaZ.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaCorrectionSolver.cpp
                         )
set(PHIST_PRECON_SOURCE ${PROJECT_SOURCE_DIR}/src/precon/phist_preconS.cpp
                        ${PROJECT_SOURCE_DIR}/src/precon/phist_preconD.cpp
                        ${PROJECT_SOURCE_DIR}/src/precon/phist_preconC.cpp
                        ${PROJECT_SOURCE_DIR}/src/precon/phist_preconZ.cpp
                        ${PROJECT_SOURCE_DIR}/src/precon/phist_ilu.cpp
                        )

include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/kernels)
include_directories(${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB})
include_directories(${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface)
include_directories(${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/krylov)
include_directories(${PROJECT_SOURCE_DIR}/src/jada)
include_directories(${PROJECT_SOURCE_DIR}/src/precon)
include_directories(${PROJECT_SOURCE_DIR}/drivers)
if (XSDK_ENABLE_CXX)
  include_directories(${PROJECT_BINARY_DIR}/cxx_bindings)
endif()
if(PHIST_KERNEL_LIB STREQUAL "magma")
  find_package(MAGMA REQUIRED)
  include_directories(${MAGMA_INCLUDE_DIRS})
  link_directories(${MAGMA_LIBRARY_DIRS})
  list(APPEND external_libs ${MAGMA_LIBRARIES})
  list(APPEND PHIST_PKG_CONFIG_REQUIRE magma)
endif()

if(PHIST_KERNEL_LIB STREQUAL "petsc")
  find_package(PETSC REQUIRED)
  include_directories(${PETSC_INCLUDE_DIRS})
  link_directories(${PETSC_LIBRARY_DIRS})
  list(APPEND external_libs ${PETSC_LIBRARIES})
  list(APPEND PHIST_PKG_CONFIG_REQUIRE PETSc)
endif()

if(PHIST_KERNEL_LIB STREQUAL "eigen")
  find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})
  add_definitions(-DEIGEN_NO_AUTOMATIC_RESIZING)
  add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_ZERO)
endif()

add_library(phist_tools ${PHIST_TOOLS_SOURCE})
add_library(phist_core ${PHIST_CORE_SOURCE})
add_library(phist_kernels_${PHIST_KERNEL_LIB} ${PHIST_KERNEL_SOURCE})
add_library(phist_solvers ${PHIST_KRYLOV_SOURCE} 
                          ${PHIST_JADA_SOURCE})
add_library(phist_precon ${PHIST_PRECON_SOURCE})

LIST(APPEND library_list phist_solvers)
LIST(APPEND library_list phist_precon)
LIST(APPEND library_list phist_core)
LIST(APPEND library_list phist_kernels_${PHIST_KERNEL_LIB})
LIST(APPEND library_list scamac)
LIST(APPEND library_list phist_tools)

target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB}                                                     phist_tools scamac)
target_link_libraries(phist_core                                          phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_precon                             phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_solvers               phist_precon phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
if(PHIST_USE_CRAFT)
  LIST(APPEND library_list craft)
  target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB}                                                     craft)
endif()

if(PHIST_KERNEL_LIB STREQUAL "ghost")
  find_package(GHOST 1.1.1 REQUIRED CONFIG)
elseif(PHIST_USE_GHOST)
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

#message(STATUS "GHOST libs: ${GHOST_LIBRARIES}")

if (PHIST_USE_MPACK)
  find_package (MPACK QUIET)
  if (MPACK_FOUND)
    target_link_libraries(phist_core ${MPACK_LIBRARIES})
    include_directories(${MPACK_INCLUDE_DIRS})
    link_directories(${MPACK_LIBRARY_DIRS})
    message(STATUS "mpack found, enabled high precision lapack!")
    if (MPACK_QD_FOUND)
      set(PHIST_HAVE_MPACK_QD 1)
    else()
      message(WARNING "found MPACK, but not MPACK_QD, which is the only 
       supported implementation right now")
    endif()
  endif()
endif()

# load parmetis for reordering builtin
if(PHIST_KERNEL_LIB STREQUAL "builtin")
  if (PHIST_USE_GRAPH_TPLS)
    find_package(ParMETIS)
    if(PARMETIS_FOUND)
      set(PHIST_HAVE_PARMETIS On)
      include_directories(${PARMETIS_INCLUDE_DIRS})
      link_directories(${PARMETIS_LIBRARY_DIRS})
      target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB} ${PARMETIS_LIBRARIES})
      message(STATUS "ParMETIS found, enabled matrix repartitioning!")
    else()
      message(STATUS "ParMETIS not found.")
    endif()

    find_package(KaHIP)
    if(KAHIP_FOUND)
      set(PHIST_HAVE_KAHIP On)
      include_directories(${KAHIP_INCLUDE_DIRS})
      link_directories(${KAHIP_LIBRARY_DIRS})
      target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB} ${KAHIP_LIBRARIES})
      message(STATUS "KAHIP found for matrix repartitioning!")
    elseif (PARMETIS_FOUND)
      message(STATUS "KaHIP not found, using ParMETIS for matrix repartitioning.")
    else()
      message(STATUS "Found neither KaHIP nor ParMETIS, no matrix repartitioning will be supported.")
    endif()

    find_package(ColPack)
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
endif()

if (PHIST_USE_LIKWID)
        set(LIKWID_INCLUDE_DIR "" CACHE PATH "Location of likwid.h")
        include_directories(${LIKWID_INCLUDE_DIR})
        get_filename_component(LIKWID_LIBDIR ${LIKWID_INCLUDE_DIR}/../lib REALPATH)
        find_library(LIKWID_LIBRARIES likwid ${LIKWID_LIBDIR})

        option(LIKWID_PERFMON "Activate function-level Likwid performance monitoring (to be used with likwid-perfctr -m)" ON)
        if (${LIKWID_PERFMON})
          message(STATUS "activated function-level Likwid performance monitoring")
        else()
          message(STATUS "disabled function-level Likwid performance monitoring")
        endif ()
        set(PHIST_HAVE_LIKWID 1)
      target_link_libraries(phist_tools ${LIKWID_LIBRARIES})

endif()


###########################################################################################
# check for required Trilinos libraries, note that you have to add the Trilinos lib/cmake #
# dir to your CMAKE_PREFIX_PATH so that CMake finds them automatically.                   #
###########################################################################################

foreach(PKG ${REQUIRE_TRILI_PKG})
  string(TOUPPER ${PKG} PKG_CAPS)
  find_package(${PKG} REQUIRED)
  message(STATUS "Look for required Trilinos package ${PKG}...")
  if (${${PKG}_FOUND})
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    target_link_libraries(phist_tools ${${PKG}_LIBRARIES})
    target_link_libraries(phist_tools ${${PKG}_TPL_LIBRARIES})
    set(PHIST_HAVE_${PKG_CAPS} "Enable required Trilinos package ${PKG}" ON)
    mark_as_advanced(${PKG}_DIR)
  else()
    message(FATAL_ERROR "required package ${PKG} was not found")
  endif()

endforeach()

#########################################################################################
# check for and enable optional libs from Trilinos (preconditioners, Belos/Anasazi)     #
#########################################################################################

if (NOT PHIST_USE_TRILINOS_TPLS)
  set(OPTIONAL_TRILI_PKG "")
endif()

foreach(PKG ${OPTIONAL_TRILI_PKG})
  string(TOUPPER ${PKG} PKG_CAPS)
  set(TPL_REPORT "optional Trilinos package ${PKG_CAPS}... ")
  find_package(${PKG} QUIET)

  if (${${PKG}_FOUND})
    list(APPEND TPL_REPORT "found")
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    target_link_libraries(phist_tools ${${PKG}_LIBRARIES})
    target_link_libraries(phist_tools ${${PKG}_TPL_LIBRARIES})
    mark_as_advanced(${PKG}_DIR)
    set(PHIST_HAVE_${PKG_CAPS} "Enable optional Trilinos package ${PKG}" ON)
  else()
    list(APPEND TPL_REPORT "not found")
  endif()
  message(STATUS ${TPL_REPORT})
endforeach()

if (PHIST_HAVE_ANASAZI)
  set(PHIST_ANASAZI_THEIR_ORTHO_MANAGER OFF CACHE LOGICAL "In the Anasazi interface, rely on the Trilinos orthogonalization routines rather than our own.")
  mark_as_advanced(PHIST_ANASAZI_THEIR_ORTHO_MANAGER)
endif()
if (PHIST_HAVE_BELOS)
  set(PHIST_BELOS_THEIR_ORTHO_MANAGER OFF CACHE LOGICAL "In the Belos interface, rely on the Trilinos orthogonalization routines rather than our own.")
  mark_as_advanced(PHIST_BELOS_THEIR_ORTHO_MANAGER)
endif()

if (PHIST_KERNEL_LIB STREQUAL "tpetra")
  if ("CUDA" IN_LIST Tpetra_TPL_LIST)
    # compiling tpetra/kokkos with CUDA is tricky and requires special flags
    # and a special wrapper script to be used as the MPI C++ compiler
    set(IS_TPETRA_CUDA 1)
    if (NOT CMAKE_CXX_COMPILER MATCHES "mpicxx" AND "$ENV{OMPI_CXX}" MATCHES "nvcc_wrapper" AND "$ENV{OMPI_MPICXX_CXXFLAGS}" MATCHES "-expt-extended-lambda")
      message(WARNING "NOTE: in order to compile Tpetra/Kokkos code with CUDA we recommend using OpenMPI"
                      "      and setting the environment variables CXX=mpicxx, OMPI_CXX=nvcc_wrapper, and OMPI_MPICXX_CXXFLAGS='-std=c++11 -expt-extended-lambda'."
                      "      The nvcc_wrapper script is installed along with Trilinos. You can unset the OMPI variables after 'make libs' to avoid compiling everytingg with nvcc.")
    endif()
  endif()
endif()


set(INCLUDE_FLAGS "")
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  set(INCLUDE_FLAGS "${INCLUDE_FLAGS} -I${dir}")
endforeach()

if(NOT GHOST_FOUND) #use MPI libs as found by GHOST in this case. otherwise, a clash between thread-safe and normal MPI will occur
  target_link_libraries(phist_tools ${MPI_CXX_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})
endif()
target_link_libraries(phist_tools m) #link the C-math library as it may not be automatically included in pure C programs

set(GCOVR_EXCLUDES)# -e src/kernels/kernels_noimpl.c -e src/tools/tools.c)

# Custom properties
set(INTEGRATION_BUILD Off CACHE BOOL "Indicates an integration build")
find_program(GCOVR_SCRIPT gcovr)
mark_as_advanced(GCOVR_SCRIPT)
find_program(CPPCHECK_EXECUTABLE cppcheck)
mark_as_advanced(CPPCHECK_EXECUTABLE)
if (${INTEGRATION_BUILD})
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} --coverage")
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --coverage")
      set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} --coverage")
    endif()
    set(GCOVR_OPTIONS -r src ${GCOVR_EXCLUDES} --xml -o coverage.xml CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS -j 24 --force --enable=all --suppress=unusedFunction --xml --xml-version=2 -q CACHE_STRING "cppcheck options")
    set(CPPCHECK_OPTIONS_POST 2> cppcheck.xml)
    set(CTEST_OPTIONS --force-new-ctest-process --verbose CACHE STRING "ctest options")
else()
    set(GCOVR_OPTIONS -r src ${GCOVR_EXCLUDES} CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS -j 24 --force --enable=all --suppress=unusedFunction CACHE STRING "cppcheck options")
    set(CPPCHECK_OPTIONS_POST)
    set(CTEST_OPTIONS --force-new-ctest-process --output-on-failure CACHE STRING "ctest options")
endif()

mark_as_advanced(CPPCHECK_OPTIONS)
mark_as_advanced(CPPCHECK_OPTIONS_POST)
mark_as_advanced(GCOVR_OPTIONS)
mark_as_advanced(CTEST_OPTIONS)

# Define includes
include_directories(lib/googletest-master-12-2015 src/ src/tools/ src/kernels/)

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
    lib/googletest-master-12-2015/gtest/gtest-all.cc
    test/mainTest.cpp
    test/tools/MatrixIO.cpp
    test/kernels/TestWithType
    test/kernels/KernelTest.cpp
    test/kernels/KernelTestWithSparseMat.cpp)
    
file(GLOB TestMatricesMM ${CMAKE_CURRENT_SOURCE_DIR}/test/matrices/*.mm)
file(GLOB TestMatricesBIN ${CMAKE_CURRENT_SOURCE_DIR}/test/matrices/*.bin)

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
# phist libraries may depend on it and the *user* cannot know
# what libraries to link. (melven)
# Additionally, this is required to make the python interface work...
target_link_libraries(phist_tools ${external_libs})

# type generic driver routines are compiled from .c or .cpp files
# which we get from .in files by preprocessing
file(GLOB TG_C_Drivers ${CMAKE_CURRENT_SOURCE_DIR}/drivers/*.c.in )
file(GLOB TG_CXX_Drivers ${CMAKE_CURRENT_SOURCE_DIR}/drivers/*.cpp.in )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/drivers)

foreach (infile ${TG_C_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    string(TOLOWER ${tc} tcl)
    set(PHIST_TG_HFILE \"phist_gen_${tcl}.h\")
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.c)
    add_executable(phist_${tc}${filebase} ${tc}${filebase}.c)
    set_property(SOURCE ${tc}${filebase}.c APPEND PROPERTY OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/phist_gen_${tcl}.h)
    target_link_libraries(phist_${tc}${filebase} ${library_list})
    install(TARGETS phist_${tc}${filebase} RUNTIME DESTINATION bin OPTIONAL)
    LIST(APPEND driver_list phist_${tc}${filebase})
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    string(TOLOWER ${tc} tcl)
    set(PHIST_TG_HFILE \"phist_gen_${tcl}.h\")
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(phist_${tc}${filebase} ${tc}${filebase}.cpp)
    set_property(SOURCE ${tc}${filebase}.cpp APPEND PROPERTY OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/phist_gen_${tcl}.h)
    target_link_libraries(phist_${tc}${filebase} ${library_list})
    install(TARGETS phist_${tc}${filebase} RUNTIME DESTINATION bin OPTIONAL)
    LIST(APPEND driver_list phist_${tc}${filebase})
  endforeach()
endforeach()

# example drivers using matfuncs

add_executable(phist_carp_cg ${PROJECT_SOURCE_DIR}/drivers/carp_cg.c)
target_link_libraries(phist_carp_cg ${library_list})

install(TARGETS phist_carp_cg RUNTIME DESTINATION bin OPTIONAL)

LIST(APPEND driver_list phist_carp_cg)

######################################################################################
# setup pkg-config and cmake package files to find phist and dependencies            #
######################################################################################

get_target_property(library_list phist_Dsubspacejada INTERFACE_LINK_LIBRARIES)

set(PHIST_PKG_CONFIG_LIBRARIES "")
foreach(lib ${library_list})
  if (IS_ABSOLUTE ${lib})
    set(PHIST_PKG_CONFIG_LIBRARIES "${PHIST_PKG_CONFIG_LIBRARIES} ${lib}")
  else()
    set(PHIST_PKG_CONFIG_LIBRARIES "${PHIST_PKG_CONFIG_LIBRARIES} -l${lib}")
  endif()
endforeach()

###########################################################################################
# workaround for tpetra and pkg-config (probably anyone using tpetra will use CMake       
# directly, but it should still work with pkg-config anyway)                              
if (PHIST_KERNEL_LIB_TPETRA)
  set(PHIST_PKG_CONFIG_LIBRARIES "${PHIST_PKG_CONFIG_LIBRARIES} -lteuchoscore")
endif()

configure_file ("${PROJECT_SOURCE_DIR}/phist.pc.in" "${PROJECT_BINARY_DIR}/phist.pc" @ONLY)
configure_file ("${PROJECT_SOURCE_DIR}/src/phist_config.h.in" "${PROJECT_BINARY_DIR}/phist_config.h")
configure_file ("${PROJECT_SOURCE_DIR}/doc/phist.doxyfile.in" "${PROJECT_BINARY_DIR}/phist.doxyfile")
include(CMakePackageConfigHelpers)
set(INCLUDE_INSTALL_DIR include)
set(LIBRARY_INSTALL_DIR lib)

configure_package_config_file("${PROJECT_SOURCE_DIR}/phistConfig.cmake.in" "${PROJECT_BINARY_DIR}/phistConfig.cmake" INSTALL_DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist
  PATH_VARS INCLUDE_INSTALL_DIR)# LIBRARY_INSTALL_DIR)
write_basic_package_version_file("${PROJECT_BINARY_DIR}/phistConfigVersion.cmake" VERSION ${PHIST_VERSION} COMPATIBILITY SameMajorVersion)

#######################################################################################
# INSTALL                                                                             #
#######################################################################################
if (${BUILD_SHARED_LIBS})
  install(TARGETS phist_kernels_${PHIST_KERNEL_LIB} EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_tools EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_core EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_precon EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_solvers EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  if(PHIST_USE_CRAFT)
    install(TARGETS craft EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  endif()

else()
  install(TARGETS phist_kernels_${PHIST_KERNEL_LIB} EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_tools EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_core EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_precon EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_solvers EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  if(PHIST_USE_CRAFT)
     install(TARGETS craft EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  endif()

endif()
install(FILES ${PROJECT_BINARY_DIR}/phistConfig.cmake        DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist)
#export(EXPORT phistLibraries FILE "${PROJECT_BINARY_DIR}/phistLibraries.cmake")
install(EXPORT phistLibraries FILE phistLibraries.cmake DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist)
install(FILES ${PROJECT_BINARY_DIR}/phistConfigVersion.cmake DESTINATION ${LIBRARY_INSTALL_DIR}/cmake/phist)
install(FILES ${PROJECT_BINARY_DIR}/phist.pc DESTINATION ${LIBRARY_INSTALL_DIR}/pkgconfig)
install(FILES ${PROJECT_BINARY_DIR}/phist_config.h DESTINATION ${INCLUDE_INSTALL_DIR})

file(GLOB PHIST_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h 
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/phist_typedefs.h 
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/*.hpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/drivers/phist_driver_utils_decl.h
                        )
install(FILES ${PHIST_HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR})

file(GLOB PYTHON_TESTS ${PROJECT_SOURCE_DIR}/python/test_*.py)
file(GLOB PYTHON_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/python/*.py)
install(FILES ${PYTHON_HEADERS} DESTINATION ${LIBRARY_INSTALL_DIR}/python)

####################
# uninstall target #
####################

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)
add_custom_target(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake") 


# test drivers
if( "${CMAKE_BUILD_TYPE}" MATCHES "Release" )
  set(PHIST_TESTS task
                  kernels
                  core
                  krylov
                  jada
                  bench
                  precon
                  tools)
else()
  set(PHIST_TESTS task
                  kernels
                  core
                  krylov
                  jada
                  precon
                  tools)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

foreach(label ${PHIST_TESTS})
  string(TOUPPER ${label} label_uc)

  file(GLOB ${label_uc}_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test/${label}/*.cpp)
  list(APPEND ${label_uc}_TEST_SOURCES ${TEST_SOURCES})
  set(TEST_EXE ${PRODUCT_NAME}-${label}-test-${CMAKE_BUILD_TYPE})
  add_executable(${TEST_EXE} ${${label_uc}_TEST_SOURCES} ${TEST_DEPENDS})
  if(PHIST_HAVE_MPI)
    set_target_properties(${TEST_EXE} PROPERTIES COMPILE_DEFINITIONS "GTEST_HAS_MPI=1;GTEST_HAS_PTHREADS=1;GTEST_HAS_EXCEPTIONS=1")
  else()
    set_target_properties(${TEST_EXE} PROPERTIES COMPILE_DEFINITIONS "GTEST_HAS_MPI=0;GTEST_HAS_PTHREADS=1;GTEST_HAS_EXCEPTIONS=1")
  endif()
  # compile tests without optimization
  set_target_properties(${TEST_EXE} PROPERTIES COMPILE_FLAGS "-O0")
  target_link_libraries(${TEST_EXE} ${library_list})
  target_link_libraries(${TEST_EXE} pthread)
  target_include_directories(${TEST_EXE} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test/${label})
  LIST(APPEND test_list ${TEST_EXE})

  # add dynmic include dependencies (which cannot be automatically detected by cmake!)
  foreach(testSrc ${${label_uc}_TEST_SOURCES})
    string(REGEX REPLACE "(With[UnAa]*lignedViews|WithB)*.cpp$" "_def.hpp" testDef ${testSrc})
    if(EXISTS ${testDef})
#      message(STATUS "Adding dynamic dependency ${testDef} to ${testSrc}")
      set_property(SOURCE ${testSrc} APPEND PROPERTY OBJECT_DEPENDS ${testDef})
    endif()
  endforeach()

endforeach()

# some example drivers which haven't made their way into the drivers/ directory
file(GLOB C_Examples ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.c )
file(GLOB CXX_Examples ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.cpp )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/examples)

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
file(GLOB TG_C_Examples ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.c.in )
file(GLOB TG_CXX_Examples ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.cpp.in )

foreach (infile ${TG_C_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    string(TOLOWER ${tc} tcl)
    set(PHIST_TG_HFILE \"phist_gen_${tcl}.h\")
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.c)
    add_executable(${tc}${filebase} ${tc}${filebase}.c)
    set_property(SOURCE ${tc}${filebase}.c APPEND PROPERTY OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/phist_gen_${tcl}.h)
    target_link_libraries(${tc}${filebase} ${library_list})
    LIST(APPEND test_list ${tc}${filebase})
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    string(TOLOWER ${tc} tcl)
    set(PHIST_TG_HFILE \"phist_gen_${tcl}.h\")
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(${tc}${filebase} ${tc}${filebase}.cpp)
    set_property(SOURCE ${tc}${filebase}.cpp APPEND PROPERTY OBJECT_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/tools/phist_gen_${tcl}.h)
    target_link_libraries(${tc}${filebase} ${library_list})
    LIST(APPEND test_list ${tc}${filebase})
  endforeach()
endforeach()

#######################################################
# TESTING
#######################################################

# Activate testing
enable_testing()

if (IS_TPETRA_CUDA)
  # Tpetra will run only on the GPU, with 12 MPI processes this gets extremely slow
  set(PHIST_MPI_TEST_NPROCS 1 2)
elseif (${PHIST_ENABLE_MPI})
  set(PHIST_MPI_TEST_NPROCS 1 2 3 12)
else()
  set(PHIST_MPI_TEST_NPROCS 1)
endif()

set(ALL_PHIST_TESTS "")

# environment variables that should be set on every MPI process
set(TEST_ENV env ASAN_OPTIONS=\"exitcode=0\" env LSAN_OPTIONS=\"exitcode=0\" env OMP_SCHEDULE=static env OMP_PROC_BIND=false env OMPI_MCA_hwloc_base_binding_policy=none)

if (PHIST_KERNEL_LIB_TPETRA)
  set(TEST_ENV ${TEST_ENV} env CUDA_LAUNCH_BLOCKING=1)
endif()
  string (REPLACE ";" " " TEST_ENV_STRING "${TEST_ENV}")
foreach(np ${PHIST_MPI_TEST_NPROCS})
  math(EXPR TEST_OMP_NUM_THREADS "12/${np}")
  set (TEST_ENV_NTHREADS ${TEST_ENV} env PHIST_NUM_THREADS=${TEST_OMP_NUM_THREADS} env OMP_NUM_THREADS=${TEST_OMP_NUM_THREADS})
  foreach(label ${PHIST_TESTS})
    set(TEST_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/test_output_np${np}/${label}")
    if( NOT PHIST_HAVE_MPI )
      set(TEST_COMMAND ${TEST_ENV_NTHREADS} 
                       ./${PRODUCT_NAME}-${label}-test-${CMAKE_BUILD_TYPE} --gtest_death_test_style=threadsafe --gtest_output=xml:${TEST_OUTPUT_DIR}/xunit.xml)
    else()
      set(TEST_COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${TEST_ENV_NTHREADS}
                       ./${PRODUCT_NAME}-${label}-test-${CMAKE_BUILD_TYPE} --gtest_death_test_style=threadsafe --gtest_output=xml:${TEST_OUTPUT_DIR}/xunit.xml)
    endif()
    set(CURRENT_TESTNAME test-${label}_np${np})
    list(APPEND ALL_PHIST_TESTS ${CURRENT_TESTNAME})
    add_test(NAME ${CURRENT_TESTNAME} COMMAND ${TEST_COMMAND})
    set_tests_properties(${CURRENT_TESTNAME} PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*PASSED[. \n]*")
    set_tests_properties(${CURRENT_TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION "[. \n]*FAILED[. \n]*")
    set_tests_properties(${CURRENT_TESTNAME} PROPERTIES TIMEOUT 1000)
    if(np EQUAL 1 AND INTEGRATION_BUILD AND PHIST_KERNEL_LIB STREQUAL "ghost" AND CUDA_LIBRARIES)
      set_tests_properties(${CURRENT_TESTNAME} PROPERTIES ENVIRONMENT GHOST_TYPE=cuda)
    else()
      set_tests_properties(${CURRENT_TESTNAME} PROPERTIES ENVIRONMENT "${MPIEXEC_ENVIRONMENT}")
    endif()
    if ( INTEGRATION_BUILD )
      # create the default test output that will appear on Jenkins if the test executable crashes
      file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
      configure_file(${PROJECT_SOURCE_DIR}/test/xunit.xml.in ${TEST_OUTPUT_DIR}/xunit.xml)
    endif()
  endforeach()
  # the task tests should be very fast (no heavy computation) - but may deadlock in case of errors!
  # JT increasing timeout because 10s is sometimes too little in debug mode (starting the executable with all the
  # dynamic libs already takes quite long, and then there's e.g. the RNG self test)
  set_tests_properties(test-task_np${np} PROPERTIES TIMEOUT 120)
endforeach()

# disable the following python tests if the GCC sanitizer is used and the kernel lib is GHOST
# because that does not work with GCC 6.1 (see issue #187)
if (PHIST_KERNEL_LIB_GHOST AND GCC_SANITIZE )
  message(WARNING "Disabling Python tests with GHOST+GCC sanitizer, cf. issue #165")
  set(DISABLE_PYTHON_TESTS 1)
endif()
# Test python interface. This requires pytest, which can be found in the python/ subdirectory.
if(NOT DISABLE_PYTHON_TESTS)
  
  set(PYSTDIR "${PROJECT_SOURCE_DIR}/python/")
  set(PYTESTDIR "${PROJECT_SOURCE_DIR}/python/pytest-3.2.2/:${PROJECT_SOURCE_DIR}/python/py-1.5.3/")
  file(COPY ${PYSTDIR}/py.test DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  
foreach(pyt ${PYTHON_TESTS})
    get_filename_component(pyt_name ${pyt} NAME_WE)
    set(pyt_name "${pyt_name}_py")
    set(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/test_output_np1/${pyt_name})
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
    set(TEST_COMMAND ${TEST_ENV} ${PYTHON_EXECUTABLE} ./py.test -s -v --junitxml ${TEST_OUTPUT_DIR}/xunit.xml --junit-prefix=${CMAKE_BUILD_TYPE} ${pyt})
    list(APPEND ALL_PHIST_TESTS ${pyt_name})
    add_test(NAME ${pyt_name} COMMAND ${TEST_COMMAND})
    set_tests_properties(${pyt_name} PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH};PYTHONPATH=${PYSTDIR}:${PYTESTDIR}:$ENV{PYTHONPATH};${MPIEXEC_ENVIRONMENT};PHIST_NUM_THREADS=12;OMP_NUM_THREADS=12")    
    if(INTEGRATION_BUILD AND PHIST_KERNEL_LIB STREQUAL "ghost" AND CUDA_LIBRARIES)
      set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "OMP_SCHEDULE=dynamic,8; GHOST_TYPE=cuda")
    endif()
  
    if( GCC_SANITIZE )
      # we need to LD_PRELOAD the respective lib<.>san.so library
      string(SUBSTRING ${GCC_SANITIZE} 0 1 GCC_SANITIZE_FIRST_CHAR)
      set(GCC_SANITIZE_LIB "lib${GCC_SANITIZE_FIRST_CHAR}san.so")
      # work around cmake bug (specifying multiple set_tests_properties(...ENVIRONMENT...)
      set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "LD_PRELOAD=${GCC_SANITIZE_LIB}:$ENV{LD_PRELOAD}")
      set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "LSAN_OPTIONS=exitcode=0 log_path=${pyt_name}_asan_log")
      set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "ASAN_OPTIONS=exitcode=0")
    endif()
  endforeach()
  if (INTEGRATION_BUILD AND CMAKE_BUILD_TYPE MATCHES "Release")
    message(STATUS "add integration tests")
    # python script that runs drivers with MPI as integration tests
    if(PHIST_HAVE_MPI)
      set(PHIST_TEST_DRIVERS_NUM_PROCS 2)
    else()
      set(PHIST_TEST_DRIVERS_NUM_PROCS 1)
    endif()
    set(PHIST_TEST_DRIVERS_NUM_THREADS 6)
    configure_file(${PROJECT_SOURCE_DIR}/python/test_phist_drivers.py.in ${CMAKE_CURRENT_BINARY_DIR}/test_phist_drivers.py)
    set(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/test_output_np${PHIST_TEST_DRIVERS_NUM_PROCS}/test_phist_drivers)
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
    set(TEST_COMMAND ${TEST_ENV} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/py.test --junitxml ${TEST_OUTPUT_DIR}/xunit.xml ./test_phist_drivers.py)
    set(pyt_name "phist_test_drivers_np${PHIST_TEST_DRIVERS_NUM_PROCS}") 
    list(APPEND ALL_PHIST_TESTS ${pyt_name})
    add_test(NAME ${pyt_name} COMMAND ${TEST_COMMAND})
    set_tests_properties(${pyt_name} PROPERTIES ENVIRONMENT "${MPIEXEC_ENVIRONMENT}; PHIST_NUM_THREADS=2; OMP_NUM_THREADS=2 OMP_SCHEDULE=dynamic,8")

    set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH}")
    set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${PYSTDIR}:${PYTESTDIR}:$ENV{PYTHONPATH}")
  endif()
endif()

# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS})

# check the installation with an example project
add_custom_target(test_install COMMAND ${PROJECT_SOURCE_DIR}/buildScripts/test_install.sh ${CMAKE_INSTALL_PREFIX})

#######################################################################################
# C++ BINDINGS                                                                        #
#######################################################################################

if (XSDK_ENABLE_CXX)
  add_subdirectory(cxx_bindings)
  # add the C++ headers as a dependency of phist_solvers so that the 'install' target works also with 'make libs'
  add_dependencies(phist_kernels_${PHIST_KERNEL_LIB} cxx_headers)
endif()

#######################################################################################
# Fortran 2003 BINDINGS                                                               #
#######################################################################################

if (XSDK_ENABLE_Fortran)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/fortran_bindings)
  add_subdirectory(fortran_bindings)
  # add the F'03 modules as a dependency of phist_solvers so that the 'install' target works also with 'make libs'
  add_dependencies(phist_kernels_${PHIST_KERNEL_LIB} fortran_modules)
endif()

#######################################################################
# DOCUMENTATION AND OTHER EXTRA TARGETS
#######################################################################

# Custom targets: doc, cppcheck, coverage, audit
# Adding optinal target "doc" to generate Doxygen documentation
find_package(Doxygen 1.7.0 QUIET)
if (DOXYGEN_FOUND)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/phist.doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    if (XSDK_ENABLE_CXX)
      add_dependencies(doc cxx_headers)
    endif()
endif()

# Add optional target "coverage" to genereate code coverage report (console / XML)
if (GCOVR_SCRIPT)
    add_custom_target(coverage ${GCOVR_SCRIPT} ${GCOVR_OPTIONS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
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
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-${CMAKE_C_COMPILER_ID}-${CMAKE_C_COMPILER_VERSION}")
set(CPACK_GENERATOR TGZ)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt" DESTINATION share/doc/phist)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/README.md"   DESTINATION share/doc/phist)
if (DOXYGEN_FOUND)
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/doc/html" DESTINATION share/doc/phist OPTIONAL)
endif()
include(CPack)

# Custom targets
add_custom_target(libs DEPENDS ${library_list})
add_custom_target(drivers DEPENDS ${driver_list})
add_custom_target(tests DEPENDS ${test_list})

