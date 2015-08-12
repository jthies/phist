# Module that checks whether MPACK is available.
# 
# Accepts the following variables:
#
# MPACK_ROOT or MPaCK_DIR: Prefix where MPACK is installed (e.g. ${HOME}/local/}.
# MLAPACK_LIB_NAME: Name of the MPACK library (default: mlapack_dd).
# MPACK_MP_LIB_NAME: Name of the library providing the basic operations (default: qd)
# MLAPACK_LIBRARY: Full path of the MPACK library
# MPACK_MP_LIBRARY: Full path of the library providing the basic operations

# Sets the following variables:
#
# MLAPACK_LIBRARY: Full path of the mlapack library.
# MPACK_FOUND: True if MPACK was found.
# MPACK_MP_LIBRARY: Full path of the multi-precision library used (e.g. gmp or qd).
# MPACK_QD_FOUND: True if basic operations are provided by QD library and it was found
# MPACK_LIBRARIES: List of all libraries needed for linking with MPACK,
# MPACK_INCLUDE_DIR
# MPACK_MP_INCLUDE_DIR
# MPACK_INCLUDE_DIRS: list of all include directories needed
# 
# Provides the following macros:
#
# find_package(MPACK)

find_path(MPACK_INCLUDE_DIR 
          NAMES mpack/mpack_config.h
          PATHS ${MPACK_DIR} ${MPACK_ROOT}
          PATH_SUFFIXES include
          NO_DEFAULT_PATH
          DOC "Include directory of MPACK")
find_path(MPACK_INCLUDE_DIR mpack/mpack_config.h)

find_path(MPACK_QD_INCLUDE_DIR qd/qd_config.h
          PATHS ${MPACK_DIR} ${MPACK_ROOT}
          PATH_SUFFIXES include
          NO_DEFAULT_PATH
          DOC "Include directory of QD")
find_path(MPACK_QD_INCLUDE_DIR qd/qd_config.h)

#check MPACK headers
include(CMakePushCheckState)
cmake_push_check_state() # Save variables
set(CMAKE_REQUIRED_INCLUDES ${CMKE_REQUIRED_INCLUDES}
                            ${MPACK_INCLUDE_DIR}
                            ${MPACK_QD_INCLUDE_DIR})
include(CheckIncludeFileCXX)
check_include_file_cxx(mpack/mpack_config.h MPACK_FOUND)
include(CheckIncludeFile)
check_include_file(qd/qd_config.h MPACK_QD_FOUND)

set(MPACK_INCLUDE_PATH ${CMAKE_REQUIRED_INCLUDES})

if (MPACK_QD_FOUND)
  set(MPACK_MP_INCLUDE_DIR ${MPACK_QD_INCLUDE_DIR})
endif()

#############
# libraries #
############# 

if (MPACK_QD_FOUND)
  set(MLAPACK_LIB_NAME mlapack_dd
    CACHE STRING "Name of the MPACK lapack library (default: mlapack_dd).")
  set(MBLAS_LIB_NAME mblas_dd
    CACHE STRING "Name of the MPACK blas library (default: mblas_dd).")
  set(MPACK_MP_LIB_NAME qd
    CACHE STRING "Name of the multi-precision library to be used with mlapack (default: qd).")
endif()

set(MLAPACK_LIBRARY MLAPACK_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the MPACK lapack library")

set(MBLAS_LIBRARY MBLAS_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the MPACK blas library")

set(MPACK_MP_LIBRARY MPACK_MP_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the MPACK library")



if(MPACK_FOUND AND MPACK_QD_FOUND)

  # search mlapack_dd library
  find_library(MLAPACK_LIBRARY ${MLAPACK_LIB_NAME}
               PATHS ${MPACK_DIR} ${MPACK_ROOT}
               PATH_SUFFIXES lib lib64
               NO_DEFAULT_PATH)
  find_library(MLAPACK_LIBRARY ${MLAPACK_LIB_NAME})

  # search mblas_dd library
  find_library(MBLAS_LIBRARY ${MBLAS_LIB_NAME}
               PATHS ${MPACK_DIR} ${MPACK_ROOT}
               PATH_SUFFIXES lib lib64
               NO_DEFAULT_PATH)
  find_library(MBLAS_LIBRARY ${MBLAS_LIB_NAME})


  # search qd library
  find_library(MPACK_MP_LIBRARY ${MPACK_MP_LIB_NAME}
               PATHS ${MPACK_DIR} ${MPACK_ROOT}
               PATH_SUFFIXES lib lib64
               NO_DEFAULT_PATH)
   find_library(MPACK_MP_LIBRARY ${MPACK_MP_LIB_NAME})
   
  # check MPACK library
#  if(MLAPACK_LIBRARY AND MPACK_MP_LIBRARY)
#    list(APPEND CMAKE_REQUIRED_LIBRARIES ${MLAPACK_LIBRARY} ${MPACK_MP_LIBRARY})
# functions in mlapack have C++ linkage, don't know if something like this would work:
#    include(CheckFunctionExists)
#    check_function_exists(Rsyevd HAVE_MPACK)
#     set(HAVE_MPACK 1)
#  endif()
  set(HAVE_MPACK 1)
endif(MPACK_FOUND AND MPACK_QD_FOUND)

if ("${MPACK_MP_LIB_NAME}" STREQUAL "qd")
  set(HAVE_MPACK_QD 1)
endif()

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "MPACK"
  DEFAULT_MSG
  MPACK_INCLUDE_DIR
  MLAPACK_LIBRARY
  MBLAS_LIBRARY
  HAVE_MPACK
)
find_package_handle_standard_args(
  "MPACK_QD"
  DEFAULT_MSG
  MPACK_QD_INCLUDE_DIR
  MPACK_MP_LIBRARY
  HAVE_MPACK_QD
)

mark_as_advanced(MPACK_INCLUDE_DIR MPACK_QD_INCLUDE_DIR MLAPACK_LIBRARY MLAPACK_LIB_NAME)
mark_as_advanced(MPACK_MP_INCLUDE_DIR MPACK_MP_LIBRARY MPACK_MP_LIB_NAME)
mark_as_advanced(MBLAS_LIB_NAME MBLAS_LIBRARY)

#restore old values
cmake_pop_check_state()

if(MPACK_FOUND)
  set(MPACK_INCLUDE_DIRS ${MPACK_INCLUDE_DIR} ${MPACK_MP_INCLUDE_DIR} CACHE STRING "required include paths for MPACK")
  set(MPACK_LIBRARIES "${MLAPACK_LIBRARY};${MBLAS_LIBRARY};${MPACK_MP_LIBRARY}" 
      CACHE FILEPATH "MPACK libraries needed for linking")
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of MPACK succeeded:\n"
    "Include directories: ${MPACK_INCLUDE_DIRS}\n"
    "Libraries: ${MPACK_LIBRARIES}\n\n")
endif(MPACK_FOUND)

