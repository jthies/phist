cmake_minimum_required (VERSION 2.6)
project (APP-SCR-TEST C Fortran CXX)

find_package(MPI REQUIRED)
if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
	link_libraries(${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()

set (CPAFT "/home/hpc/unrz/unrza354/Projects/ESSEX/GHOST/CPLIB/CpAftlib/")
include_directories(${CPAFT})
link_directories("${CPAFT}")
set(LIBS ${LIBS} ${CPAFT})
link_libraries (cpaft)

add_executable(minimal.bin minimal.cpp )


