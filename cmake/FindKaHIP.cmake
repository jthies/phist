# Module that checks whether KaHIP is available (the Karlsruhe partitioning package).
# 
# Accepts the following variables:
#
# KAHIP_ROOT: Prefix where KaHIP is installed.
# KAHIP_LIB_NAME: Name of the KaHIP library (default: parmetis).
# KAHIP_LIBRARY: Full path of the KaHIP library

# Sets the following variables:
#
# KAHIP_LIBRARY: Full path of the KaHIP library.
# KAHIP_FOUND: True if KaHIP was found.
# KAHIP_LIBRARIES: List of all libraries needed for linking with KaHIP,
# 
# Provides the following macros:
#
# find_package(KaHIP)

find_path(KAHIP_INCLUDE_DIR kaHIP_interface.h
          PATHS ${KAHIP_DIR} ${KAHIP_ROOT}
          PATH_SUFFIXES include kahip KaHIP kaHIP
          NO_DEFAULT_PATH
          DOC "Include directory of KaHIP")
find_path(KAHIP_INCLUDE_DIR kaHIP_interface.h
          PATH_SUFFIXES include kahip KaHIP kaHIP)

set(KAHIP_LIB_NAME kahip
    CACHE STRING "Name of the KaHIP library (default: kahip).")
set(KAHIP_LIBRARY KaHIP_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the KaHIP library")


#check KaHIP headers
include(CMakePushCheckState)
cmake_push_check_state() # Save variables
include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${KAHIP_INCLUDE_DIR})
check_include_file_cxx(kaHIP_interface.h KAHIP_FOUND)

if(KAHIP_FOUND)
  set(KaHIP_INCLUDE_PATH ${CMAKE_REQUIRED_INCLUDES})

  # search KaHIP library
  find_library(KAHIP_LIBRARY ${KAHIP_LIB_NAME}
               PATHS ${KAHIP_DIR} ${KAHIP_ROOT}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH)
  find_library(KAHIP_LIBRARY ${KAHIP_LIB_NAME})

  # check KaHIP library
  if(KAHIP_LIBRARY)
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${KAHIP_LIBRARY})
# cmake 3.8 hangs on this one, and since it's not declared extern "C" in KaHIP,
# we can't use check_function_exists as we do for ParMETIS
#    include(CheckCXXSymbolExists)
#    check_cxx_symbol_exists( kaffpaE ${KAHIP_LIBRARY} HAVE_KAHIP)
    set(HAVE_KAHIP ON)
  endif(KAHIP_LIBRARY)
endif(KAHIP_FOUND)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "KaHIP"
  DEFAULT_MSG
  KAHIP_INCLUDE_DIR
  KAHIP_LIBRARY
  HAVE_KAHIP
)

mark_as_advanced(KAHIP_INCLUDE_DIR KAHIP_LIBRARY KAHIP_LIB_NAME)

#restore old values
cmake_pop_check_state()


if(KAHIP_FOUND)
  set(KAHIP_INCLUDE_DIRS ${KAHIP_INCLUDE_DIR})
  set(KAHIP_LIBRARIES "${KAHIP_LIBRARY};${METIS_LIBRARY}" 
      CACHE FILEPATH "KaHIP libraries needed for linking")
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of KaHIP succeded:\n"
    "Include directory: ${KAHIP_INCLUDE_DIRS}\n"
    "Library directory: ${KAHIP_LIBRARIES}\n\n")
endif(KAHIP_FOUND)
