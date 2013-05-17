##############################################################################
#
# CMake variable for use by Trilinos/Tpetra clients. 
#
# Do not edit: This file was generated automatically by CMake.
#
##############################################################################

#
# Make sure CMAKE_CURRENT_LIST_DIR is usable
#

IF (NOT DEFINED CMAKE_CURRENT_LIST_DIR)
  GET_FILENAME_COMPONENT(_THIS_SCRIPT_PATH ${CMAKE_CURRENT_LIST_FILE} PATH)
  SET(CMAKE_CURRENT_LIST_DIR ${_THIS_SCRIPT_PATH})
ENDIF()


## ---------------------------------------------------------------------------
## Compilers used by Trilinos/Tpetra build
## ---------------------------------------------------------------------------

SET(Tpetra_CXX_COMPILER "/usr/bin/mpiCC")

SET(Tpetra_C_COMPILER "/usr/bin/mpicc")

SET(Tpetra_FORTRAN_COMPILER "/usr/bin/mpif90")


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Tpetra build
## ---------------------------------------------------------------------------

## Set compiler flags, including those determined by build type
SET(Tpetra_CXX_FLAGS "  -g -fPIC -O3")

SET(Tpetra_C_FLAGS "  -g -fPIC -O3")

SET(Tpetra_FORTRAN_FLAGS " -O3")

## Extra link flags (e.g., specification of fortran libraries)
SET(Tpetra_EXTRA_LD_FLAGS "-g -fPIC")

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty. 
SET(Tpetra_SHARED_LIB_RPATH_COMMAND "")
SET(Tpetra_BUILD_SHARED_LIBS "OFF")

SET(Tpetra_LINKER /usr/bin/ld)
SET(Tpetra_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths 
## ---------------------------------------------------------------------------

## List of package include dirs
SET(Tpetra_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../include")

## List of package library paths
SET(Tpetra_LIBRARY_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../lib")

## List of package libraries
SET(Tpetra_LIBRARIES "tpetraext;tpetrainout;tpetra;epetra;kokkosdisttsqr;kokkosnodetsqr;kokkoslinalg;kokkosnodeapi;kokkos;tpi;teuchosremainder;teuchosnumerics;teuchoscomm;teuchosparameterlist;teuchoscore;teuchosremainder;teuchosnumerics;teuchoscomm;teuchosparameterlist;teuchoscore")

## Specification of directories for TPL headers
SET(Tpetra_TPL_INCLUDE_DIRS "/usr/include")

## Specification of directories for TPL libraries
SET(Tpetra_TPL_LIBRARY_DIRS "")

## List of required TPLs
SET(Tpetra_TPL_LIBRARIES "/usr/lib/liblapack.so;/usr/lib/libblas.so;/usr/lib/libblas.so;/usr/lib/i386-linux-gnu/libpthread.so")

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

SET(Tpetra_MPI_LIBRARIES "")
SET(Tpetra_MPI_LIBRARY_DIRS "")
SET(Tpetra_MPI_INCLUDE_DIRS "")
SET(Tpetra_MPI_EXEC "MPI_EXEC-NOTFOUND")
SET(Tpetra_MPI_EXEC_MAX_NUMPROCS "")
SET(Tpetra_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables 
## ---------------------------------------------------------------------------

## The packages enabled for this project
SET(Tpetra_PACKAGE_LIST "Tpetra;Epetra;KokkosClassic;ThreadPool;Teuchos;TeuchosRemainder;TeuchosNumerics;TeuchosComm;TeuchosParameterList;TeuchosCore")

## The TPLs enabled for this project
SET(Tpetra_TPL_LIST "LAPACK;BLAS;MPI;Pthread")


# Import Trilinos targets
IF(NOT Trilinos_TARGETS_IMPORTED)
  SET(Trilinos_TARGETS_IMPORTED 1)
  INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../Trilinos/TrilinosTargets.cmake")
ENDIF()

