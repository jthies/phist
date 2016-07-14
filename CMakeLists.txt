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

set (SCR "/home/hpc/unrz/unrza354/local/scr-modif-ulfm/")
include_directories ("${SCR}/include/")
link_directories ("${SCR}/lib/")
set(LIBS ${LIBS} ${SCR}/lib/libscr.so)

link_libraries (scr)
add_definitions(-DSCR)
#find_library(SCR)
#if(SCR_FOUND)
#	message (STATUS "SCR pacage is found ${SCR_INCLUDE_DIRS} ")
#	include_directories(${SCR_INCLUDE_DIRS})
#	set(LIBS ${LIBS} ${SCR_LIBRARIES})
#else()
#	message (STATUS "SCR package is NOT found")
#endif()



message (STATUS "LIBS ${LIBS}")
message (STATUS "INCLUDE ${INCLUDE_DIRECTORIES}")

add_library(cpaft SHARED Checkpoint.cpp cp_options.cpp aft.cpp aft_macros.h )

