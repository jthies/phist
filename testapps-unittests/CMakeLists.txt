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

if(GHOST_TEST)
	find_package(GHOST REQUIRED)
	if(GHOST_FOUND)
		message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
		include_directories(${GHOST_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${GHOST_LIBRARIES})
		link_libraries(${GHOST_LIBRARIES})
	else()
		message (STATUS "GHOST package is NOT found")
	endif()
endif()

find_package(CRAFT REQUIRED)
if(CRAFT_FOUND)
	message (STATUS "CRAFT package is found ${CRAFT_INCLUDE_PATH}")
	message (STATUS "CRAFT_FLAGS: ${CRAFT_FLAGS}")
	include_directories(${CRAFT_INCLUDE_PATH})
	set(LIBS ${LIBS} ${CRAFT_LIBRARIES})
	link_libraries(${CRAFT_LIBRARIES})
	set( CMAKE_CXX_FLAGS ${CRAFT_FLAGS})
else()
	message (STATUS "CRAFT package is NOT found")
endif()





if(GHOST_TEST)

		find_package(ESSEX-PHYSICS REQUIRED)
		if (ESSEX-PHYSICS_FOUND)
    		message(STATUS "Found ESSEX-Physics ${ESSEX-PHYSICS_VERSION}")
				set(LIBS ${LIBS} ${ESSEX-PHYSICS_LIBRARIES})
				include_directories(${ESSEX-PHYSICS_INCLUDE_DIRS})
				set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GHOST_C_FLAGS}")
				set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GHOST_CXX_FLAGS}")
				set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GHOST_Fortran_FLAGS}")
				link_directories("${ESSEX-PHYSICS_LIBRARIES}")
				link_libraries(${ESSEX-PHYSICS_LIBRARIES})
    		add_definitions(-DHAVE_ESSEX_PHYSICS)
		else()
    		message(STATUS "Could not find ESSEX-Physics! Some examples will be skipped!")
		endif()

	find_path(ESSEXAMPLES_INC REQURIED)
	find_path(ESSEXAMPLES_LIB REQURIED)
	message (STATUS "ESSEXAMPLES_LIB is found ${ESSEXAMPLES_LIB}")
	include_directories(${ESSEXAMPLES_INC})
	link_directories("${ESSEXAMPLES_LIB}")
	link_libraries (essexamples)
endif()

add_subdirectory(minimal)
add_subdirectory(arrayTest)
if(GHOST_TEST)
	add_subdirectory(ghostTest)
endif()

