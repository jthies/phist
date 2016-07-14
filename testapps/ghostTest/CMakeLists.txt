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

find_package(GHOST REQUIRED)
    if(GHOST_FOUND)
		message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
		include_directories(${GHOST_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${GHOST_LIBRARIES})
		link_libraries(${GHOST_LIBRARIES})
	else()
		message (STATUS "GHOST package is NOT found")
endif()

find_package(ESSEX-PHYSICS REQUIRED)
	if(ESSEX-PHYSICS_FOUND)
		message (STATUS "GHOST pacage is found ${ESSEX-PHYSICS_INCLUDE_DIRS} ")
		include_directories(${ESSEX-PHYSICS_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${ESSEX-PHYSICS_LIBRARIES})
		link_libraries(${ESSEX-PHYSICS_LIBRARIES})
	else()
			message (STATUS "ESSEX-PHYSICS package is NOT found")
endif()



set (CPAFT "/home/hpc/unrz/unrza354/Projects/ESSEX/GHOST/CPLIB/CpAftlib/")
include_directories(${CPAFT})
link_directories("${CPAFT}")
set(LIBS ${LIBS} ${CPAFT})
link_libraries (cpaft)

set (GHOSTSANDBOX "/home/hpc/unrz/unrza354/Projects/ESSEX/GHOST/ghost-sandbox")
set (ESSEXAMPLES_INC "/home/hpc/unrz/unrza354/Projects/ESSEX/GHOST/ghost-sandbox/ghost-apps/common")
set (ESSEXAMPLES_LIB "/home/hpc/unrz/unrza354/Projects/ESSEX/GHOST/ghost-sandbox/build-ulfm-aft/ghost-apps/common/")
include_directories(${ESSEXAMPLES_INC})
link_directories("${ESSEXAMPLES_LIB}")
set(LIBS ${LIBS} ${ESSEXAMPLES_LIB})
link_libraries (essexamples)


add_executable(ghostTest.bin ghostTest.cpp )

