cmake_minimum_required(VERSION 2.8.10)

# some guidelines for extending/modifying the PHIST CMakeLists.txt file. 
# Things may not be handled this way everywhere yet, but you should keep these
# points in mind.
#
# - anything prefixed by PHIST_ENABLE or PHIST_USE will enter the CMake Cache but not the config file,
#   so it won't be directly visible in the source code
# - anything prefixed by PHIST_HAVE can be seen int he phist_config.h file
# - other things should *not* go into the cache unless they come from another library via find_package.
#
# - try to mark things that are very detailed (like details on the timing behavior etc.) as ADVANCED so that
#   the CMake cache seen when first using ccmake stays reasonable

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
#list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Trilinos")

# Define some custom variables
project(phist C CXX Fortran)
SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)
set(PHIST_VERSION 0.9.2)
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
option(PHIST_BUILD_MIC "Build for Intel MIC" OFF)

option(PHIST_HAVE_SSE  "enable SSE SIMD intrinsics in PHIST code" ON)
option(PHIST_HAVE_AVX  "enable AVX SIMD intrinsics in PHIST code" OFF)
option(PHIST_HAVE_AVX2 "enable AVX2 SIMD intrinsics in PHIST code" OFF)

# reproducible random number generator?
# TODO - this is not quite working yet, once it is, make it the default
option(PHIST_BUILTIN_RNG
        "use random number generator provided by PHIST, which may be solower than the one implemented by the kernel lib. On the other hand, it gives reproducible results across runs, independent of the exact runtime setup."
        ON)

# use GHOST if found (must be ON if ghost is the kernel lib)
option(PHIST_USE_GHOST
        "try to find the GHOST library even if it is not the kernel lib" 
        ON)

# try to find MPACK for high precision lapack routines
option(PHIST_USE_MPACK
        "try to find the MPACK library" 
        ON)

# try to find third-party solver packages (depending on kernel lib),
# for instance for epetra/tpetra and GHOST we have an interface to
# Anasazi and Belos.
option(PHIST_USE_SOLVER_TPLS
        "try to find supported third-party libraries that implement iterative methods for the kernel library used."
        ON)

# try to find third-party preconditioning packages (depending on kernel lib),
# for instance if the kernel lib is epetra, look for Ifpack and ML.
option(PHIST_USE_PRECON_TPLS
        "try to find supported third-party libraries that implement preconditioners for the kernel library used."
        ON)

# try to find third-party partitioning packages (depending on kernel lib),
# for instance if the kernel lib is epetra, look for Isorropia, builtin supports ParMETIS 
# and ColPack etc.
option(PHIST_USE_GRAPH_TPLS
        "try to find supported third-party libraries that implement graph partitioning and reordering for the kernel library used."
        ON)

# this option allows to disable all Trilinos libraries at once, overriding
# any of the choices PHIST_USE_*_TPLS
option(PHIST_USE_TRILINOS_TPLS
        "enable/disable using any optional Trilinos libraries like Zoltan, Belos, Kokkos etc."
        ON)


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
find_path(LAPACKE_INCL NAMES "lapacke.h" "mkl_lapacke.h" $ENV{MKLROOT}/include)
if(LAPACKE_INCL STREQUAL "LAPACKE_INCL-NOTFOUND")
  message(FATAL_ERROR "LAPACKe include file not found. Add the proper path to CMAKE_PREFIX_PATH.")
else()
  message(STATUS "LAPACKe include file found in: ${LAPACKE_INCL}")
  include_directories(${LAPACKE_INCL})
endif()

# Try to find an MKL header in the default location. Do this first, so
# in case the MKL is available but no Intel compiler is used, not
# finding LAPACK(e)/BLAS later on is not an error.
if (NOT MKL_INCLUDE_DIR)
  message(STATUS "Try to find MKL in default paths and via MKLROOT")
  find_path(MKL_INCLUDE_DIR mkl.h $ENV{MKLROOT}/include)
endif()
if (MKL_INCLUDE_DIR)
  set(PHIST_HAVE_MKL 1)
  include_directories(${MKL_INCLUDE_DIR})
endif()

# For non-Intel compilers, check for additional LAPACK(e)/BLAS libraries.
# If the MKL is present, not finding the required libs is not an error.
# TODO Always prefer MKL over others?
if(NOT CMAKE_C_COMPILER_ID STREQUAL "Intel")
# TODO: this assumes that if you use the intel compilers, you also have MKL (!)
#=======================================================================
  if(PHIST_HAVE_MKL)
    # TODO write a module for that
#-----------------------------------------------------------------------
    # interface
    unset(MKL_INTERFACE_LIBS CACHE)
    find_library(MKL_INTERFACE_LIBS NAMES mkl_intel_lp64 PATHS "$ENV{MKLROOT}/lib/intel64")
    if(MKL_INTERFACE_LIBS STREQUAL "MKL_INTERFACE_LIBS-NOTFOUND")
      message(FATAL_ERROR "Could not find MKL interface library. Is the MKLROOT environment variable set properly?")
    else()
      message(STATUS "MKL interface library found: ${MKL_INTERFACE_LIBS}")
      LIST(APPEND external_libs ${MKL_INTERFACE_LIBS})
    endif()
#-----------------------------------------------------------------------
    # threading or sequential
    unset(MKL_THREADING_LIBS CACHE)
    if( PHIST_ENABLE_OPENMP )
      # linking the GNU variant, we shouldn't need libiomp5.
      find_library(MKL_THREADING_LIBS NAMES mkl_gnu_thread PATHS "$ENV{MKLROOT}/lib/intel64")
    else()
      find_library(MKL_THREADING_LIBS NAMES mkl_sequential PATHS "$ENV{MKLROOT}/lib/intel64")
    endif()
    if(MKL_THREADING_LIBS STREQUAL "MKL_THREADING_LIBS-NOTFOUND")
      message(FATAL_ERROR "Could not find MKL threading library. Is the MKLROOT environment variable set properly?")
    else()
      message(STATUS "MKL threading library found: ${MKL_THREADING_LIBS}")
      LIST(APPEND external_libs ${MKL_THREADING_LIBS})
    endif()
#-----------------------------------------------------------------------
    # core libs
    unset(MKL_COMPUTING_LIBS CACHE)
    find_library(MKL_COMPUTING_LIBS NAMES mkl_core PATHS "$ENV{MKLROOT}/lib/intel64")
    if(MKL_COMPUTING_LIBS STREQUAL "MKL_COMPUTING_LIBS-NOTFOUND")
      message(FATAL_ERROR "Could not find MKL computing library. Is the MKLROOT environment variable set properly?")
    else()
      message(STATUS "MKL computing library found: ${MKL_COMPUTING_LIBS}")
      LIST(APPEND external_libs ${MKL_COMPUTING_LIBS})
    endif()
  else()
#=======================================================================
    # First search external non-system paths to allow the use of 
    # non-default libraries. Only if nothing is found there also search
    # system paths.
    # Maybe even add a special variable for LAPACK, LAPACKe and BLAS
    # paths.

    # LAPACKe
    unset(LAPACKE_LIBS CACHE)
    find_library(LAPACKE_LIBS NAMES "liblapacke.so" "liblapacke.a" NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    if(LAPACKE_LIBS STREQUAL "LAPACKE_LIBS-NOTFOUND")
      # we could add NO_CMAKE_PATH and NO_CMAKE_ENVIRONMENT_PATH because we already searched there
      find_library(LAPACKE_LIBS NAMES "liblapacke.so" "liblapacke.a")
    endif()
    if(LAPACKE_LIBS STREQUAL "LAPACKE_LIBS-NOTFOUND")
      message(FATAL_ERROR "LAPACKe libraries not found. Add the proper path to CMAKE_PREFIX_PATH.")
    else()
      message(STATUS "LAPACKe libraries found: ${LAPACKE_LIBS}")
      LIST(APPEND external_libs ${LAPACKE_LIBS})
    endif()
#-----------------------------------------------------------------------
    # LAPACK
    unset(LAPACK_LIBS CACHE)
    find_library(LAPACK_LIBS NAMES "liblapack.so" "liblapack.a" NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    if(LAPACK_LIBS STREQUAL "LAPACK_LIBS-NOTFOUND")
      find_library(LAPACK_LIBS NAMES "liblapack.so" "liblapack.a")
    endif()
    if(LAPACK_LIBS STREQUAL "LAPACK_LIBS-NOTFOUND")
      message(FATAL_ERROR "LAPACK libraries not found. Add the proper path to CMAKE_PREFIX_PATH.")
    else()
      message(STATUS "LAPACK libraries found: ${LAPACK_LIBS}")
      LIST(APPEND external_libs ${LAPACK_LIBS})
    endif()
#-----------------------------------------------------------------------
    # BLAS
    unset(BLAS_LIBS CACHE)
    find_library(BLAS_LIBS NAMES "libblas.so" "libblas.a" NO_CMAKE_SYSTEM_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    if(BLAS_LIBS STREQUAL "BLAS_LIBS-NOTFOUND")
      find_library(BLAS_LIBS NAMES "libblas.so" "libblas.a")
    endif()
    if(BLAS_LIBS STREQUAL "BLAS_LIBS-NOTFOUND")
      message(FATAL_ERROR "BLAS libraries not found. Add the proper path to CMAKE_PREFIX_PATH.")
    else()
      message(STATUS "BLAS libraries found: ${BLAS_LIBS}")
      LIST(APPEND external_libs ${BLAS_LIBS})
    endif()
#=======================================================================
  endif()
endif()
########################################################################

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

  if (CMAKE_C_COMPILER_VERSION VERSION_LESS 4.8)
    set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter --coverage")
    set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter --coverage")
    set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -DTESTING -Wall -fcheck=all -fstack-protector-all --coverage")
  else ()
    set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -fstack-protector-all -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-local-typedefs --coverage")
    set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -fstack-protector-all -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-local-typedefs --coverage")
    set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -DTESTING -Wall -fcheck=all -fstack-protector-all --coverage")
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

  if (PHIST_BUILD_MIC)
      set (CMAKE_C_FLAGS_RELEASE        "-O3 -no-prec-div -mmic")
      set (CMAKE_CXX_FLAGS_RELEASE      "-O3 -no-prec-div -mmic")
      set (CMAKE_Fortran_FLAGS_RELEASE  "-O3 -no-prec-div -mmic")
  else()
      set (CMAKE_C_FLAGS_RELEASE        "-O3 -no-prec-div -xHOST")
      set (CMAKE_CXX_FLAGS_RELEASE      "-O3 -no-prec-div -xHOST")
      set (CMAKE_Fortran_FLAGS_RELEASE  "-O3 -no-prec-div -xHOST")
  endif()

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -debug")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -debug")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -debug")

  set (CMAKE_C_FLAGS_DEBUG        "-check-pointers=rw -fstack-protector -w2 -debug -traceback -DTESTING")
  set (CMAKE_CXX_FLAGS_DEBUG      "-check-pointers=rw -fstack-protector -w2 -debug -traceback -DTESTING")
  set (CMAKE_Fortran_FLAGS_DEBUG  "-check all -debug -traceback -DTESTING")
  if (PHIST_BUILD_MIC)
      set (CMAKE_C_FLAGS_DEBUG        "-mmic -check-pointers=rw -check=uninit,stack -fstack-protector -w2 -debug -traceback -DTESTING")
      set (CMAKE_CXX_FLAGS_DEBUG      "-mmic -check-pointers=rw -check=uninit,stack -fstack-protector -w2 -debug -traceback -DTESTING")
      set (CMAKE_Fortran_FLAGS_DEBUG  "-mmic -check all -debug -traceback -DTESTING")
  endif()
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

elseif (CMAKE_C_COMPILER_ID STREQUAL "PGI")

  message(STATUS "using PGI compiler flags")

  set (CMAKE_C_FLAGS        "-c99")
  set (CMAKE_CXX_FLAGS      "-std=c++0x")
  set (CMAKE_Fortran_FLAGS  "-Mpreprocess")
  if( PHIST_ENABLE_OPENMP )
    set(PHIST_HAVE_OPENMP 1)
    set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -mp=align")
    set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -mp=align")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -mp=align")
  endif()

  if(PHIST_HAVE_MKL)
    set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS}       -DMKL_LP64")
    set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS}     -DMKL_LP64")
  endif()

  # -ffast kills high precision stuff!
  set (CMAKE_C_FLAGS_RELEASE        "-O4 -fast -Mvect=levels:10")
  set (CMAKE_CXX_FLAGS_RELEASE      "-O2")
  set (CMAKE_Fortran_FLAGS_RELEASE  "-O4 -fast -Mvect=levels:10")

  set (CMAKE_C_FLAGS_RELWITHDEBINFO       "${CMAKE_C_FLAGS_RELEASE}       -g")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO     "${CMAKE_CXX_FLAGS_RELEASE}     -g")
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -g")

  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g")#-O0 -g -DTESTING -Mbounds -Mdwarf3")
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0")#-O0 -g -DTESTING -Mbounds -Mdwarf3")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")#-O0 -g -DTESTING -Mbounds -Mdwarf3")

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
CHECK_CXX_SOURCE_COMPILES("
  int main(int argc, char** argv)
  {
    thread_local int x;
    return 0;
  }
" PHIST_HAVE_CXX11_THREADLOCAL)
CHECK_CXX_SOURCE_COMPILES("
  template<typename T>
  class C {
  public:
  C(){}
  C(C&&)=default;
  };
  int main(int argc, char** argv)
  {
    C<int> c;
  }
" PHIST_HAVE_CXX11_MOVEDEFAULT)
if(PHIST_HAVE_CXX11_LAMBDAS OR PHIST_HAVE_CXX11_THREADLOCAL)
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
  set(MPIEXEC_ENVIRONMENT "OMPI_MCA_hwloc_base_binding_policy=none" CACHE STRING "")
endif()

##############################################################################
# setup which kernel library is being used.
##############################################################################

set(PHIST_KERNEL_LIB "builtin" CACHE STRING "Select underlying kernel libraries for phist (builtin|ghost|epetra|tpetra|magma)")
message(STATUS "kernel lib " ${PHIST_KERNEL_LIB} " ...")

set(REQUIRE_TRILI_PKG,"")
set(OPTIONAL_TRILI_PKG,"")

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
  if (PHIST_USE_SOLVER_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Anasazi")
    list(APPEND OPTIONAL_TRILI_PKG "Belos")
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

  if (PHIST_USE_SOLVER_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Anasazi")
    list(APPEND OPTIONAL_TRILI_PKG "Belos")
  endif()
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
  set(PHIST_SDMATS_ROW_MAJOR 0 CACHE STRING "store sdMats in row-major order (not supported by our algorithms..." )

  set(PHIST_SELL_C -1 CACHE INTEGER "chunk size C for the SELL-C-sigma format (-1: automatic choice)")
  set(PHIST_SELL_SIGMA -1 CACHE INTEGER "sorting scope sigma for the SELL-C-sigma format (-1: automatic choice)")

  if (PHIST_USE_SOLVER_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Anasazi")
    list(APPEND OPTIONAL_TRILI_PKG "Belos")
  endif()
  if (PHIST_USE_GRAPH_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Zoltan")
  endif()
  list(APPEND OPTIONAL_TRILI_PKG "Kokkos")
  list(APPEND OPTIONAL_TRILI_PKG "Teuchos")

elseif (PHIST_KERNEL_LIB STREQUAL "builtin")
  set(PHIST_KERNEL_LIB_BUILTIN 1)
  set(PHIST_MVECS_ROW_MAJOR 1)
  option(PHIST_KERNEL_LIB_BUILTIN_PIN_THREADS "Try to pin the threads to cores using sched_setaffinity (based on simple heuristics; for Intel use env KMP_AFFINITY=\"granularity=fine,compact,1,0\"" On)
  if(PHIST_HAVE_AVX2)
    message(STATUS "Activated experimental high-precision builtin kernel routines!")
    set(PHIST_HIGH_PRECISION_KERNELS 1)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -mavx2 -mfma")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mavx2 -mfma")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx -mavx2 -mfma")

    if (NOT MPACK_FOUND)
      message(WARNING "did not find mlapack_dd library, as you are compiling with high precision kernels, you may want to provide this library."
                      "It can be obtained from http://mplapack.sourceforge.net/ and should be compiled with support for double-double (dd) core operations."
                      "If it is already installed and you want to use it, set MPACK_DIR to provide a hint where it can be found.")
    endif()
  endif()

  if( NOT PHIST_HAVE_MPI )
    message(FATAL_ERROR "Builtin kernel lib requires MPI!")
  endif()
elseif (PHIST_KERNEL_LIB STREQUAL "magma")
  set(PHIST_KERNEL_LIB_MAGMA 1)
  set(PHIST_MVECS_ROW_MAJOR 0)
  set(PHIST_SDMATS_ROW_MAJOR 0)
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set or not recognized" )
endif ()
set(PHIST_BENCH_LARGE_N 100000000 CACHE STRING "Dimension for benchmarking tests (e.g. large enough to make micro benchmarks memory-bound). Should be a multiple of 8.")

option(PHIST_TIMINGS_FULL_TRACE "Try to print full call trace in the TimeMonitor/performance measurements (Warning: may generate *a lot* of output)." Off)

# even without the essex-physics libs we include certain matfuncs for the
# test cases SpinSZ and Graphene for testing purposes.
find_package(ESSEX-PHYSICS CONFIG QUIET)
if (ESSEX-PHYSICS_FOUND)
  message(STATUS "Found ESSEX-Physics ${ESSEX-PHYSICS_VERSION}")
  set(PHIST_HAVE_ESSEX_PHYSICS 1)
  include_directories(${ESSEX-PHYSICS_INCLUDE_DIRS})
else()
  message(STATUS "ESSEX-Physics not found, will use supplementary implementation")
  include_directories(${PROJECT_SOURCE_DIR}/drivers/matfuncs)
endif()

#############################################################################
# use 32-bit integers as global indices with Epetra by default. This is     #
# because e.g. Isorropia/Zoltan and our own HYMLS require it. For larger    #
# runs with Epetra it has to be switched off before compiling PHIST.        #
#############################################################################
if( "${PHIST_KERNEL_LIB}" MATCHES "epetra" )
  set(PHIST_FORCE_INT_GIDX_DEFAULT ON)
else()
  set(PHIST_FORCE_INT_GIDX_DEFAULT OFF)
endif()

option(PHIST_FORCE_INT_GIDX "Use 32-bit global indices for compatiblity with other libs" 
${PHIST_FORCE_INT_GIDX_DEFAULT})

##############################################################################
# enable single precision and complex valued calculation if supported
##############################################################################
if( "${PHIST_KERNEL_LIB}" MATCHES "tpetra|ghost|magma" )
  set(PHIST_SUPPORT_SP ON)
  set(PHIST_SUPPORT_CMPLX ON)
else()
  set(PHIST_SUPPORT_SP OFF)
  set(PHIST_SUPPORT_CMPLX OFF)
endif()

option(PHIST_ENABLE_SP "Enable single precision functions and drivers" ${PHIST_SUPPORT_SP})
option(PHIST_ENABLE_CMPLX "Enable complex valued functions and drivers" ${PHIST_SUPPORT_CMPLX})

if (PHIST_ENABLE_SP)
  if(NOT PHIST_SUPPORT_SP)
    message(WARNING "Enabled PHIST_ENABLE_SP, but ${PHIST_KERNEL_LIB} kernel lib does not support single-precision calculations!")
  endif()
  set(PHIST_HAVE_SP 1)
endif ()

if (PHIST_ENABLE_CMPLX)
  if(NOT PHIST_SUPPORT_CMPLX)
    message(WARNING "Enabled PHIST_ENABLE_CMPLX, but ${PHIST_KERNEL_LIB} kernel lib does not support complex-valued calculations!")
  endif()
  set(PHIST_HAVE_CMPLX 1)
endif()

# these specific versions of type-generic drivers have to be created
if (PHIST_ENABLE_SP)
  if (PHIST_ENABLE_CMPLX)
    set(TG_Types S D C Z)
  else ()
    set(TG_Types S D)
  endif ()
else ()
  if (PHIST_ENABLE_CMPLX)
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
                          ${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface/phist_rcp_helpers.cpp
                          )
  if(PHIST_HIGH_PRECISION_KERNELS)
    set (PHIST_KERNEL_SOURCE ${PHIST_KERNEL_SOURCE}
                          ${PROJECT_SOURCE_DIR}/src/kernels/builtin/prec_helpers.c
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
                           ${PROJECT_SOURCE_DIR}/src/kernels/phist_mpi_kernels.cpp
                           ${PROJECT_SOURCE_DIR}/src/kernels/trilinos_interface/phist_rcp_helpers.cpp)
endif()
if (NOT ESSEX-PHYSICS_FOUND)
  LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/drivers/matfuncs/graphene.c
                                  ${PROJECT_SOURCE_DIR}/drivers/matfuncs/SpinChainSZ.c
                                  ${PROJECT_SOURCE_DIR}/drivers/matfuncs/part_tools.c)
endif()

# repartitioning with Isorropia/Zoltan
if (PHIST_KERNEL_LIB STREQUAL "epetra")
  LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/epetra/hg_repart.cpp)
endif()

# sopme `derived' kernel functions, e.g. mvec_my_length is implemented via mvec_get_map/map_get_local_lenth etc.
# and common things like benchmarks
LIST(APPEND PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/common/kernels_common_impl.cpp
                                ${PROJECT_SOURCE_DIR}/src/kernels/common/stream_bench.c)

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
                       ${PROJECT_SOURCE_DIR}/src/core/phist_kpm.cpp
                       )
set (PHIST_KRYLOV_SOURCE ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_arnoldi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_lanczos.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedgmres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_blockedminres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belos.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_anasazi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_carp_cg.cpp
                         )
set (PHIST_JADA_SOURCE   ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOpts.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jdqr.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_schur_decomp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejada.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_transform_searchspace.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaCorrectionSolver.cpp
                         )
set(PHIST_PRECON_SOURCE ${PROJECT_SOURCE_DIR}/src/precon/phist_precon.cpp
                        ${PROJECT_SOURCE_DIR}/src/precon/phist_ilu.cpp
                        )
set (PHIST_SUPPORT_SOURCE  ${PROJECT_SOURCE_DIR}/src/support/phist_feastCorrectionSolver.cpp
                           ${PROJECT_SOURCE_DIR}/src/support/phist_iter_op.cpp
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
include_directories(${PROJECT_SOURCE_DIR}/src/support)
include_directories(${PROJECT_SOURCE_DIR}/drivers)

add_library(phist_tools ${PHIST_TOOLS_SOURCE})
add_library(phist_core ${PHIST_CORE_SOURCE})
add_library(phist_kernels_${PHIST_KERNEL_LIB} ${PHIST_KERNEL_SOURCE})
add_library(phist_solvers ${PHIST_KRYLOV_SOURCE} 
                          ${PHIST_JADA_SOURCE})
add_library(phist_precon ${PHIST_PRECON_SOURCE})
add_library(phist_support ${PHIST_SUPPORT_SOURCE})

LIST(APPEND library_list phist_support)
LIST(APPEND library_list phist_solvers)
LIST(APPEND library_list phist_precon)
LIST(APPEND library_list phist_core)
LIST(APPEND library_list phist_kernels_${PHIST_KERNEL_LIB})
LIST(APPEND library_list phist_tools)

target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB}                                                     phist_tools)
target_link_libraries(phist_core                                          phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_precon                             phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_solvers               phist_precon phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)
target_link_libraries(phist_support phist_solvers phist_precon phist_core phist_kernels_${PHIST_KERNEL_LIB} phist_tools)

if (MPACK_FOUND)
      target_link_libraries(phist_core ${MPACK_LIBRARIES})
endif()
if (ESSEX-PHYSICS_FOUND)
  target_link_libraries(phist_kernels_${PHIST_KERNEL_LIB} ${ESSEX-PHYSICS_LIBRARIES})
endif()

if(PHIST_KERNEL_LIB STREQUAL "ghost")
  find_package(GHOST REQUIRED CONFIG)
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
  find_package (MPACK)
  if (MPACK_FOUND)
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

if(PHIST_KERNEL_LIB STREQUAL "magma")
  find_package(MAGMA REQUIRED)
  include_directories(${MAGMA_INCLUDE_DIRS})
  link_directories(${MAGMA_LIBRARY_DIRS})
  list(APPEND external_libs ${MAGMA_LIBRARIES})
  list(APPEND PHIST_PKG_CONFIG_REQUIRE magma)
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
      message(STATUS "ParMETIS found, enabled matrix reordering!")
    else()
      message(STATUS "ParMETIS not found, no matrix reordering!")
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

###########################################################################################
# check for required Trilinos libraries, note that you have to add the Trilinos lib/cmake #
# dir to your CMAKE_PREFIX_PATH so that CMake finds them automatically.                   #
###########################################################################################

foreach(PKG ${REQUIRE_TRILI_PKG})
  string(TOUPPER ${PKG} PKG_CAPS)
  message(STATUS "look for ${PKG} (required)")
  find_package(${PKG} REQUIRED)
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
  set(OPTIONAL_TRILINOS_PKG "")
endif()

foreach(PKG ${OPTIONAL_TRILI_PKG})
  string(TOUPPER ${PKG} PKG_CAPS)
  message(STATUS "look for ${PKG} (optional)")
  find_package(${PKG} QUIET)
  if (${${PKG}_FOUND})
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    target_link_libraries(phist_tools ${${PKG}_LIBRARIES})
    target_link_libraries(phist_tools ${${PKG}_TPL_LIBRARIES})
    mark_as_advanced(${PKG}_DIR)
    set(PHIST_HAVE_${PKG_CAPS} "Enable optional Trilinos package ${PKG}" ON)
  endif()

endforeach()


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
  option(PHIST_PERFCHECK_SEPARATE_OUTPUT "write perfcheck output to separate files for each MPI rank (useful for analysing heterogenous runs)" OFF)
  option(PHIST_PERFCHECK_REALISTIC "Switch from ideal performance predictions to more realistic ones (e.g. for views)!" OFF)
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
    # added
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
  install(TARGETS phist_precon EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_solvers EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_support EXPORT phistLibraries LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
else()
  install(TARGETS phist_kernels_${PHIST_KERNEL_LIB} EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_tools EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_core EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_precon EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_solvers EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
  install(TARGETS phist_support EXPORT phistLibraries ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR})
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

  file(GLOB ${label_uc}_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test/${label}/*.cpp)
  list(APPEND ${label_uc}_TEST_SOURCES ${TEST_SOURCES})
  set(TEST_EXE ${PRODUCT_NAME}-${label}-test-${CMAKE_BUILD_TYPE})
  add_executable(${TEST_EXE} ${${label_uc}_TEST_SOURCES} ${TEST_DEPENDS})
  set_target_properties(${TEST_EXE} PROPERTIES COMPILE_DEFINITIONS "GTEST_HAS_MPI=1;GTEST_HAS_PTHREADS=1;GTEST_HAS_EXCEPTIONS=1")
  # compile tests without optimization
  set_target_properties(${TEST_EXE} PROPERTIES COMPILE_FLAGS "-O0")
  target_link_libraries(${TEST_EXE} ${library_list})
  target_link_libraries(${TEST_EXE} pthread)
  target_include_directories(${TEST_EXE} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test/${label})
  LIST(APPEND test_list ${TEST_EXE})

  # add dynmic include dependencies (which cannot be automatically detected by cmake!)
  foreach(testSrc ${${label_uc}_TEST_SOURCES})
    string(REGEX REPLACE "(With[UnAa]*lignedViews)*.cpp$" "_def.hpp" testDef ${testSrc})
    if(EXISTS ${testDef})
#      message(STATUS "Adding dynamic dependency ${testDef} to ${testSrc}")
      set_property(SOURCE ${testSrc} APPEND PROPERTY OBJECT_DEPENDS ${testDef})
    endif()
  endforeach()

endforeach()

# some example drivers which haven't made their way into the drivers/ directory
file(GLOB C_Examples ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.c )
file(GLOB CXX_Examples ${CMAKE_CURRENT_SOURCE_DIR}/examples/*/*.cpp )

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

#

# Activate testing
enable_testing()

if (${PHIST_ENABLE_MPI})
  set(PHIST_MPI_TEST_NPROCS 1 2 3 12)
else()
  set(PHIST_MPI_TEST_NPROCS 1)
endif()

foreach(np ${PHIST_MPI_TEST_NPROCS})
  math(EXPR TEST_OMP_NUM_THREADS "12/${np}")
  foreach(label ${PHIST_TESTS})
    set(TEST_COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} env LSAN_OPTIONS=\"exitcode=0\" env PHIST_NUM_THREADS=${TEST_OMP_NUM_THREADS} env OMP_NUM_THREADS=${TEST_OMP_NUM_THREADS} env OMP_SCHEDULE=static
    ./${PRODUCT_NAME}-${label}-test-${CMAKE_BUILD_TYPE} --gtest_death_test_style=threadsafe --gtest_output=xml:test_output_np${np}/${label}/xunit.xml)
    add_test(NAME test-${label}_np${np} COMMAND ${TEST_COMMAND})
    set_tests_properties(test-${label}_np${np} PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*PASSED[. \n]*")
    set_tests_properties(test-${label}_np${np} PROPERTIES FAIL_REGULAR_EXPRESSION "[. \n]*FAILED[. \n]*")
    set_tests_properties(test-${label}_np${np} PROPERTIES TIMEOUT 500)
    set_tests_properties(test-${label}_np${np} PROPERTIES ENVIRONMENT ${MPIEXEC_ENVIRONMENT})
  endforeach()
  # the task tests should be very fast (no heavy computation) - but may deadlock in case of errors!
  # JT increasing timeout because 10s is sometimes too little in debug mode (starting the executable with all the
  # dynamic libs already takes quite long, and then there's e.g. the RNG self test)
  set_tests_properties(test-task_np${np} PROPERTIES TIMEOUT 120)
endforeach()

# disable the following python tests if the GCC sanitizer is used and the kernel lib is Tpetra
# because that does not work (cf. #165)
#if (PHIST_KERNEL_LIB_TPETRA AND GCC_SANITIZE )
#  message(WARNING "Disabling Python tests with Tpetra+GCC sanitizer, cf. issue #165")
#  set(DISABLE_PYTHON_TESTS 1)
#endif()
# Test python interface if this is an integration build.
# The reason why this test is only included in integration builds is that it requires
# py.test, which is included in phist but currently with a hard-coded interpreter from
# the module system (anaconda module). We can't just load the module because it brings its own
# MPI installation which breaks the rest of PHIST.
if(INTEGRATION_BUILD)
  configure_file(${PROJECT_SOURCE_DIR}/python/py.test ${CMAKE_CURRENT_BINARY_DIR}/ COPYONLY)
  # note: we can do system-specific stuff like chmod here because we know that it's an integration build
  # and no one should run it except for ourselves
  execute_process(COMMAND "chmod" "u+x" "${CMAKE_CURRENT_BINARY_DIR}/py.test")
  foreach(pyt ${PYTHON_TESTS})
    get_filename_component(pyt_name ${pyt} NAME_WE)
    set(pyt_name "${pyt_name}_py")
    message(STATUS "add python test '${pyt_name}'")
    set(TEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/test_output_np1/${pyt_name})
    file(MAKE_DIRECTORY ${TEST_OUTPUT_DIR})
    set(TEST_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/py.test --junitxml ${TEST_OUTPUT_DIR}/xunit.xml ${pyt})
    add_test(NAME ${pyt_name} COMMAND ${TEST_COMMAND})
    set_tests_properties(${pyt_name} PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{LD_LIBRARY_PATH} PYTHONPATH=${PYSTDIR}:$ENV{PYTHONPATH} PHIST_NUM_THREADS=12 OMP_NUM_THREADS=12 OMP_SCHEDULE=dynamic,8")
  
    if( GCC_SANITIZE )
      # we need to LD_PRELOAD the respective lib<.>san.so library
      string(SUBSTRING ${GCC_SANITIZE} 0 1 GCC_SANITIZE_FIRST_CHAR)
      set(GCC_SANITIZE_LIB "lib${GCC_SANITIZE_FIRST_CHAR}san.so")
      # work around cmake bug (specifying multiple set_tests_properties(...ENVIRONMENT...)
      set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "LD_PRELOAD=${GCC_SANITIZE_LIB}:$ENV{LD_PRELOAD}")
      set_property(TEST ${pyt_name} APPEND PROPERTY ENVIRONMENT "LSAN_OPTIONS=\"exitcode=0\"")
    endif()
  endforeach()
endif()

# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} ${CTEST_OPTIONS})

# Custom targets: doc, cppcheck, coverage, audit
# Adding optinal target "doc" to generate Doxygen documentation
find_package(Doxygen 1.7.0)
if (DOXYGEN_FOUND)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/phist.doxyfile WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
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

