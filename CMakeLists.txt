cmake_minimum_required (VERSION 2.6)
project (CPAFTLIB)
find_package(MPI REQUIRED)

if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()

find_package(GHOST REQUIRED)
if(GHOST_FOUND)
	message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
	include_directories(${GHOST_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${GHOST_LIBRARIES})
else()
	message (STATUS "GHOST package is NOT found")
endif()


find_path(PHIST REQUIRED)
include_directories ("${PHIST}/include/")
link_directories ("${PHIST}/lib/")


add_library(cpaft SHARED Checkpoint.cpp cp_options.cpp aft_macros.h CpPhistMvec.h CpPhistMvec.cpp CpPhistMvec_decl.h CpPhistMvec_def.hpp)
