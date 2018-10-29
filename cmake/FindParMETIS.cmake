# Module that checks whether ParMETIS is available.
#
# Accepts the following variables:
#
# TPL_ParMETIS_DIR: Prefix where ParMETIS is installed.
# TPL_ParMETIS_INCLUDE_DIRS: list of directories of installed ParMETIS
# TPL_ParMETIS_LIBRARIES: list of libraries of installed ParMETIS.
# METIS_LIB_NAME: Name of the METIS library (default: metis).
# PARMETIS_LIB_NAME: Name of the ParMETIS library (default: parmetis).
# METIS_LIBRARY: Full path of the METIS library.
# PARMETIS_LIBRARY: Full path of the ParMETIS library

# Sets the following variables:
#
# PARMETIS_INCLUDE_DIRS
# METIS_LIBRARY: Full path of the METIS library (unless the user gave us TPL_ParMETIS_LIBRARIES).
# PARMETIS_LIBRARY: Full path of the ParMETIS library ( " ).
# PARMETIS_FOUND: True if ParMETIS was found.
# PARMETIS_LIBRARIES: List of all libraries needed for linking with ParMETIS,
#
# Provides the following macros:
#
# find_package(ParMETIS)

if (TPL_ParMETIS_INCLUDE_DIRS)
  find_path(PARMETIS_INCLUDE_DIR parmetis.h
          PATHS ${TPL_ParMETIS_INCLUDE_DIRS}
          NO_DEFAULT_PATH
          DOC "Include directory of ParMETIS"
          REQUIRED
          )
elseif(TPL_ParMETIS_DIR)
  find_path(PARMETIS_INCLUDE_DIR parmetis.h
          PATHS ${TPL_ParMETIS_DIR}
          PATH_SUFFIXES include parmetis include/parmetis
          NO_DEFAULT_PATH
          REQUIRED)
else()
  find_path(PARMETIS_INCLUDE_DIR parmetis.h
          PATH_SUFFIXES include parmetis)

endif()
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})

set(METIS_LIB_NAME metis
    CACHE STRING "Name of the METIS library (default: metis).")
set(PARMETIS_LIB_NAME parmetis
    CACHE STRING "Name of the ParMETIS library (default: parmetis).")
set(METIS_LIBRARY METIS_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the METIS library")
set(PARMETIS_LIBRARY ParMETIS_LIBRARY-NOTFOUND
    CACHE FILEPATH "Full path of the ParMETIS library")


#check METIS and ParMETIS headers
include(CMakePushCheckState)
cmake_push_check_state() # Save variables
include(CheckIncludeFile)
set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${MPI_INCLUDE_PATH} ${PARMETIS_INCLUDE_DIR})
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${MPI_DUNE_COMPILE_FLAGS}")
check_include_file(metis.h METIS_FOUND)
check_include_file(parmetis.h PARMETIS_FOUND)

if(PARMETIS_FOUND)
  set(ParMETIS_INCLUDE_PATH ${CMAKE_REQUIRED_INCLUDES})
  set(ParMETIS_COMPILE_FLAGS "${CMAKE_REQUIRED_FLAGS} -DENABLE_PARMETIS=1")

  if (NOT TPL_ParMETIS_LIBRARIES)
    # search METIS library
    find_library(METIS_LIBRARY metis
               PATHS ${TPL_ParMETIS_DIR} ${PARMETIS_DIR} ${PARMETIS_ROOT}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH)
    find_library(METIS_LIBRARY metis)

    # search ParMETIS library
    find_library(PARMETIS_LIBRARY parmetis
               PATHS ${TPL_ParMETIS_DIR} ${PARMETIS_DIR} ${PARMETIS_ROOT}
               PATH_SUFFIXES lib
               NO_DEFAULT_PATH)
    find_library(PARMETIS_LIBRARY parmetis)
    set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY} CACHE FILEPATH "ParMETIS library list.")
  else()
    set(PARMETIS_LIBRARIES ${TPL_ParMETIS_LIBRARIES} CACHE FILEPATH "ParMETIS library list.")
  endif()
  # check ParMETIS library
  if(PARMETIS_LIBRARIES)
    list(APPEND CMAKE_REQUIRED_LIBRARIES ${PARMETIS_LIBRARIES})
    include(CheckFunctionExists)
    check_function_exists(parmetis_v3_partkway HAVE_PARMETIS)
  endif(PARMETIS_LIBRARIES)
endif(PARMETIS_FOUND)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "ParMETIS"
  DEFAULT_MSG
  PARMETIS_INCLUDE_DIRS
  PARMETIS_LIBRARIES
  HAVE_PARMETIS
)

mark_as_advanced(PARMETIS_INCLUDE_DIR METIS_LIBRARY PARMETIS_LIBRARY METIS_LIB_NAME PARMETIS_LIB_NAME)

#restore old values
cmake_pop_check_state()


if(PARMETIS_FOUND)
  if (NOT PARMETIS_LIBRARIES)
    set(PARMETIS_LIBRARIES "${PARMETIS_LIBRARY};${METIS_LIBRARY}" 
        CACHE FILEPATH "ParMETIS libraries needed for linking")
  endif()
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ParMETIS succeded:\n"
    "Include directory: ${PARMETIS_INCLUDE_DIR}\n"
    "Libraries: ${PARMETIS_LIBRARIES}\n\n")
endif(PARMETIS_FOUND)
