cmake_minimum_required (VERSION 2.6)
project (GHOST-TEST C Fortran CXX)


###### OPTIONS ######
option(GHOST_TEST
				"In order to compile GHOST_TEST, ESSEXAMPELS-library will be required"
		OFF)
######################


find_package(MPI REQUIRED)
if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
	link_libraries(${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()

find_package(GHOST REQUIRED)
    if(GHOST_FOUND)
		message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
		include_directories(${GHOST_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${GHOST_LIBRARIES})
		link_libraries(${GHOST_LIBRARIES})
	else()
		message (STATUS "GHOST package is NOT found")
endif()

find_package(CPAFT REQUIRED)
if(CPAFT_FOUND)
	message (STATUS "CPAFT package is found ${CPAFT_INCLUDE_PATH}")
	message (STATUS "CPAFT_FLAGS: ${CPAFT_FLAGS}")
	include_directories(${CPAFT_INCLUDE_PATH})
	set(LIBS ${LIBS} ${CPAFT_LIBRARIES})
	link_libraries(${CPAFT_LIBRARIES})
	set( CMAKE_CXX_FLAGS ${CPAFT_FLAGS})
else()
	message (STATUS "CPAFT package is NOT found")
endif()

if(GHOST_TEST)
	find_path(ESSEXAMPLES_INC REQURIED)
	find_path(ESSEXAMPLES_LIB REQURIED)
	include_directories(${ESSEXAMPLES_INC})
	link_directories("${ESSEXAMPLES_LIB}")
	link_libraries (essexamples)
endif()

add_subdirectory(minimal)
add_subdirectory(arrayTest)
if(GHOST_TEST)
	add_subdirectory(ghostTest)
endif()

