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

if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()


if(GHOST_CP)
	find_package(GHOST REQUIRED)
	if(GHOST_FOUND)
		message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
		include_directories(${GHOST_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${GHOST_LIBRARIES})
		add_definitions(-DGHOST_CP)
	else()
		message (STATUS "GHOST package is NOT found")
	endif()
endif()


if(PHIST_CP)
	find_path(PHIST REQUIRED)
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
	add_definitions(-DPHIST_CP)
endif()

if(SCR)
	find_package(MPI REQUIRED)
	if(SCR_FOUND)
		message (STATUS "SCR package is found ${SCR_INCLUDE_PATH}")
		include_directories(${SCR_INCLUDE_PATH})
		set(LIBS ${LIBS} ${SCR_LIBRARIES})
	else()
		message (STATUS "SCR package is NOT found")
	endif()
endif()

set ( CPAFT_SOURCE 	${SRC_DIR}/Checkpoint.cpp
		${SRC_DIR}/cp_options.cpp)

if(AFT)
	set ( CPAFT_SOURCE 	${CPAFT_SOURCE}
						${SRC_DIR}/aft.cpp
						${SRC_DIR}/aft_macros.h)
endif()	



add_library(cpaft SHARED ${CPAFT_SOURCE} ${CPPHISTMVEC_SRC} ${CPPHISTSDMAT_SRC})



