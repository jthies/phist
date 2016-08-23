cmake_minimum_required (VERSION 2.6)
project (CPAFTLIB)
set (SRC_DIR ${CMAKE_SOURCE_DIR})


###### OPTIONS ######
option(PHIST_CP
		"enables/disables the checkpoints of PHIST data types"
		OFF)
option(GHOST_CP
		"enables/disables the checkpoints of PHIST data types"
		OFF)
option(SCR
		"enables/disables the usage of SCR-library for node-level C/R"
		OFF)
option(AFT
		"Automatic Fault Tolerance: This feature can only be anabled in case of ULFM-MPI"
		OFF)


###### FIND LIBRARIES ######

find_package(MPI REQUIRED)
if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()


if(GHOST_CP)
	find_package(GHOST REQUIRED)
	add_definitions(-DGHOST_CP)
	if(GHOST_FOUND)
		message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
		include_directories(${GHOST_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${GHOST_LIBRARIES})
	else()
		message (STATUS "GHOST package is NOT found")
	endif()
endif()


if(PHIST_CP)
	find_path(PHIST REQUIRED)
	add_definitions(-DPHIST_CP)
	include_directories ("${PHIST}/include/")
	link_directories ("${PHIST}/lib/")


	set (CPPHISTMVEC_DIR 	${SRC_DIR}/cpTypes/cpPhistMvec/)
	include_directories ("${CPPHISTMVEC_DIR}")
	set (CPPHISTMVEC_SRC 	${CPPHISTMVEC_DIR}/CpPhistMvec.h
							${CPPHISTMVEC_DIR}/CpPhistMvec_decl.h
							${CPPHISTMVEC_DIR}/CpPhistMvec.cpp
							${CPPHISTMVEC_DIR}/CpPhistMvec_def.hpp)

	set (CPPHISTSDMAT_DIR 	${SRC_DIR}/cpTypes/cpPhistSdMat/)
	include_directories ("${CPPHISTSDMAT_DIR}")
	set (CPPHISTSDMAT_SRC 	${CPPHISTSDMAT_DIR}/CpPhistSdMat.h
							${CPPHISTSDMAT_DIR}/CpPhistSdMat_decl.h
							${CPPHISTSDMAT_DIR}/CpPhistSdMat.cpp
							${CPPHISTSDMAT_DIR}/CpPhistSdMat_def.hpp)
endif()


if(SCR)
	find_library(SCR_PATH NAMES scr)

	if(SCR_PATH-NOTFOUND)
		message (STATUS "SCR package is NOT found")
	else()
		message (STATUS "SCR lib path ${SCR_PATH}")
		set ( SCR_LIBRARY_PATH "${SCR_PATH}/lib/")
		message (STATUS "library path ${SCR_LIBRARY_PATH}")
		set ( SCR_INCLUDE_PATH "${SCR_PATH}/include/" )
		message (STATUS "include path ${SCR_INCLUDE_PATH}")
		include_directories("${SCR_INCLUDE_PATH}")
		link_directories ("${SCR_LIBRARY_PATH}")
		set(LIBS ${LIBS} ${SCR_LIBRARY_PATH})
		add_definitions(-DSCR)
	endif()
endif()

set (CPAFT_SOURCE 	${SRC_DIR}/Checkpoint.cpp
		${SRC_DIR}/cp_options.cpp
		)

if(AFT)
	find_file(MPI_EXT_PATH "mpi-ext.h" HINTS ${MPI_INCLUDE_PATH})
	message (STATUS "MPI_EXT_PATH: ${MPI_EXT_PATH}")
	if("${MPI_EXT_PATH}" STREQUAL "MPI_EXT_PATH-NOTFOUND")
		message (FATAL_ERROR "mpi-ext.h is not found. Compile with ULFM-MPI for using AFT feature.")
	endif()
		set ( CPAFT_SOURCE 	${CPAFT_SOURCE}
							${SRC_DIR}/aft.cpp
							${SRC_DIR}/aft_macros.h)
endif()	

add_library(cpaft SHARED ${CPAFT_SOURCE} ${CPPHISTMVEC_SRC} ${CPPHISTSDMAT_SRC})
set(LIBRARY_INSTALL_DIR lib)

configure_file("${SRC_DIR}/cpaft-config.cmake.in" "${SRC_DIR}/cpaft-config.cmake" @ONLY)


##### INSTALL TARGET #####
install ( TARGETS cpaft EXPORT cpaftLibrary LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR} )
install ( FILES ${SRC_DIR}/cpaft-config.cmake DESTINATION ${LIBRARY_INSTALL_DIR} )

file(GLOB CPAFT_HEADERS ${CMAKE_SOURCE_DIR}/*.h
												${CMAKE_SOURCE_DIR}/*.hpp	
												)

file(GLOB CPAFT_TYPE_SRC 
					${CMAKE_SOURCE_DIR}/cpTypes	
		)	
set(INCLUDE_INSTALL_DIR include)
install ( FILES ${CPAFT_HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR} )
install ( DIRECTORY ${CPAFT_TYPE_SRC} DESTINATION ${INCLUDE_INSTALL_DIR} )

##### UNINSTALL TARGET #####
configure_file(
				    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
						 "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	    				IMMEDIATE @ONLY)

add_custom_target(uninstall
		COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)


