set(LIKWID_HOME "$ENV{LIKWID_HOME}" CACHE STRING "base directory of Likwid installation")
if(LIKWID_HOME)
  message(STATUS "Searching for Likwid in ${LIKWID_HOME}")
  find_library(LIKWID_LIBRARIES likwid PATHS ${LIKWID_HOME} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  if(LIKWID_LIBRARIES)
    set(LIKWID_Found true)
    set(LIKWID_INCLUDE_DIRS ${LIKWID_HOME}/include)
    set(LIKWID_LIBRARY_DIRS ${LIKWID_HOME}/lib)
  endif()
endif()
