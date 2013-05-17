##############################################################################
#
# CMake variable for use by Trilinos/Epetra clients. 
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
## Compilers used by Trilinos/Epetra build
## ---------------------------------------------------------------------------

SET(Epetra_CXX_COMPILER "/usr/bin/mpiCC")

SET(Epetra_C_COMPILER "/usr/bin/mpicc")

SET(Epetra_FORTRAN_COMPILER "/usr/bin/mpif90")


## ---------------------------------------------------------------------------
## Compiler flags used by Trilinos/Epetra build
## ---------------------------------------------------------------------------

## Set compiler flags, including those determined by build type
SET(Epetra_CXX_FLAGS "  -g -fPIC -O3")

SET(Epetra_C_FLAGS "  -g -fPIC -O3")

SET(Epetra_FORTRAN_FLAGS " -O3")

## Extra link flags (e.g., specification of fortran libraries)
SET(Epetra_EXTRA_LD_FLAGS "-g -fPIC")

## This is the command-line entry used for setting rpaths. In a build
## with static libraries it will be empty. 
SET(Epetra_SHARED_LIB_RPATH_COMMAND "")
SET(Epetra_BUILD_SHARED_LIBS "OFF")

SET(Epetra_LINKER /usr/bin/ld)
SET(Epetra_AR /usr/bin/ar)

## ---------------------------------------------------------------------------
## Set library specifications and paths 
## ---------------------------------------------------------------------------

## List of package include dirs
SET(Epetra_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../include")

## List of package library paths
SET(Epetra_LIBRARY_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../lib")

## List of package libraries
SET(Epetra_LIBRARIES "epetra")

## Specification of directories for TPL headers
SET(Epetra_TPL_INCLUDE_DIRS "")

## Specification of directories for TPL libraries
SET(Epetra_TPL_LIBRARY_DIRS "")

## List of required TPLs
SET(Epetra_TPL_LIBRARIES "/usr/lib/liblapack.so;/usr/lib/libblas.so;/usr/lib/libblas.so")

## ---------------------------------------------------------------------------
## MPI specific variables
##   These variables are provided to make it easier to get the mpi libraries
##   and includes on systems that do not use the mpi wrappers for compiling
## ---------------------------------------------------------------------------

SET(Epetra_MPI_LIBRARIES "")
SET(Epetra_MPI_LIBRARY_DIRS "")
SET(Epetra_MPI_INCLUDE_DIRS "")
SET(Epetra_MPI_EXEC "MPI_EXEC-NOTFOUND")
SET(Epetra_MPI_EXEC_MAX_NUMPROCS "")
SET(Epetra_MPI_EXEC_NUMPROCS_FLAG "")

## ---------------------------------------------------------------------------
## Set useful general variables 
## ---------------------------------------------------------------------------

## The packages enabled for this project
SET(Epetra_PACKAGE_LIST "Epetra")

## The TPLs enabled for this project
SET(Epetra_TPL_LIST "LAPACK;BLAS;MPI")


# Import Trilinos targets
IF(NOT Trilinos_TARGETS_IMPORTED)
  SET(Trilinos_TARGETS_IMPORTED 1)
  INCLUDE("${CMAKE_CURRENT_LIST_DIR}/../Trilinos/TrilinosTargets.cmake")
ENDIF()

