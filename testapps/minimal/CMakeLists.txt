cmake_minimum_required (VERSION 2.6)
project (APP-SCR-TEST C Fortran CXX)
set (SRC_DIR ${CMAKE_SOURCE_DIR})

find_package(MPI REQUIRED)
if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
	link_libraries(${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()

find_path(CPAFT REQIRED)
include_directories ("${CPAFT}/include")
link_directories ("${CPAFT}/lib")
link_libraries (cpaft)

add_executable(minimal.bin minimal.cpp )

