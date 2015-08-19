# Module that checks whether ColPack is available.
# 
# Accepts the following variables:
#
# COLPACK_ROOT: Prefix where ColPack is installed.
# COLPACK_LIB_NAME: Name of the ColPack library (default: ColPack).
# COLPACK_LIBRARY: Full path of the ColPack library

# Sets the following variables:
#
# COLPACK_LIBRARY: Full path of the ColPack library.
# COLPACK_FOUND: True if ColPack was found.
# COLPACK_LIBRARIES: List of all libraries needed for linking with ColPack,
# 
# Provides the following macros:
#
# find_package(ColPack)

find_path(COLPACK_INCLUDE_DIR ColPackHeaders.h
          PATHS ${COLPACK_DIR} ${COLPACK_ROOT}
          PATH_SUFFIXES include ColPack
          DOC "Include directory of ColPack")

set(COLPACK_LIB_NAME ColPack
          PATHS ${COLPACK_DIR} ${COLPACK_ROOT}
          PATH_SUFFIXES include include/ColPack
          CACHE STRING "Name of the ColPack library (default: ColPack).")
set(COLPACK_LIBRARY ColPack_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the ColPack library")

#check ColPack headers
include(CMakePushCheckState)
cmake_push_check_state() # Save variables
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${COLPACK_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
check_include_file_cxx(ColPackHeaders.h COLPACK_FOUND)
if(COLPACK_FOUND)
  set(ColPack_INCLUDE_PATH ${CMAKE_REQUIRED_INCLUDES})
  set(ColPack_COMPILE_FLAGS "${CMAKE_REQUIRED_FLAGS}")

  # search ColPack library
  find_library(COLPACK_LIBRARY ColPack
               PATHS ${COLPACK_DIR} ${COLPACK_ROOT}
               PATH_SUFFIXES lib lib64
               NO_DEFAULT_PATH)
  find_library(COLPACK_LIBRARY ColPack)

  # check ColPack library
#  if(COLPACK_LIBRARY)
#    list(APPEND CMAKE_REQUIRED_LIBRARIES ${COLPACK_LIBRARY})
#    include(CheckFunctionExists)
#    check_function_exists(ColPack_v3_partkway HAVE_COLPACK)
#  endif(COLPACK_LIBRARY)
  set(HAVE_COLPACK 1)
endif(COLPACK_FOUND)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ColPack"
  DEFAULT_MSG
  COLPACK_INCLUDE_DIR
  COLPACK_LIBRARY
  HAVE_COLPACK
)

mark_as_advanced(COLPACK_INCLUDE_DIR COLPACK_LIBRARY COLPACK_LIB_NAME)

#restore old values
cmake_pop_check_state()


if(COLPACK_FOUND)
  set(COLPACK_INCLUDE_DIRS ${COLPACK_INCLUDE_DIR})
  set(COLPACK_LIBRARIES "${COLPACK_LIBRARY}" 
      CACHE FILEPATH "ColPack libraries needed for linking")
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ColPack succeded:\n"
    "Include directory: ${COLPACK_INCLUDE_DIRS}\n"
    "Library directory: ${COLPACK_LIBRARIES}\n\n")
endif(COLPACK_FOUND)
