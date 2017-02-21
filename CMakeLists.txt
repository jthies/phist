cmake_minimum_required (VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)
project (CRAFTLIB)
set (PROJECT_DIR ${CMAKE_SOURCE_DIR})
set (PROJECT_SRC_DIR ${PROJECT_DIR}/src/)
set (PROJECT_INC_DIR ${PROJECT_DIR}/include/)
set (PROJECT_TOOLS_DIR ${PROJECT_DIR}/tools/)

###### OPTIONS ######
option(PHIST_CP
		"enables/disables the checkpoints of PHIST data types"
		OFF)
option(GHOST_CP
		"enables/disables the checkpoints of PHIST data types"
		OFF)
option(MKL_CP
		"enables/disables the usage of MKL data types for C/R"
		ON)
option(SCR
		"enables/disables the usage of SCR-library for node-level C/R"
		OFF)
option(AFT
		"Automatic Fault Tolerance: This feature can only be anabled in case of ULFM-MPI"
		OFF)
option(MPIIO
		"If SCR is not defined, use MPI-IO for simple POD, POD-array checkpoints"
		ON)


set (SOURCE 	
		src/checkpoint.cpp
		src/cpOptions.cpp
		src/cpArray.cpp
		src/cpPOD.cpp
		src/cpHelperFuncs.cpp
		src/craftConf.cpp
		src/timing.cpp
		)

set (HEADERS
		include/checkpoint.hpp
		include/cpArray.hpp
		include/cpBase.hpp
		include/cpPOD.hpp
		include/cpHelperFuncs.hpp
		include/cpEnum.h
		include/cpOptions.h
		include/dataType.h
		include/craftConf.h
		include/timing.h
    include/addedCpTypes.hpp
		)

set(TOOLS
    tools/rm-cache-ckpts.sh
    tools/do-rm-cache-ckpts.sh
)


include_directories(${HEADERS} ${PROJECT_DIR} ${PROJECT_INC_DIR})

###### FIND LIBRARIES ######

find_package(MPI REQUIRED)
if(MPI_FOUND)
	message (STATUS "MPI package is found ${MPI_INCLUDE_PATH}")
	include_directories(${MPI_INCLUDE_PATH})
	set(LIBS ${LIBS} ${MPI_LIBRARIES})
else()
	message (STATUS "MPI package is NOT found")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(AFT)
	find_file(MPI_EXT_PATH "mpi-ext.h" HINTS ${MPI_INCLUDE_PATH})
	message (STATUS "MPI_EXT_PATH: ${MPI_EXT_PATH}")
	if("${MPI_EXT_PATH}" STREQUAL "MPI_EXT_PATH-NOTFOUND")
		message (FATAL_ERROR "mpi-ext.h is not found. Compile with ULFM-MPI for using AFT feature.")
	endif()
		set ( SOURCE 	${SOURCE}
							src/aft.cpp)
		set ( HEADERS 	${HEADERS}
							include/aft_macros.h)
endif()	


if(SCR)
  if(MPIIO)
   message( FATAL_ERROR "MPI_IO can not be used with SCR.")
  endif()

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

if(MPIIO)
  if(SCR)
   message( FATAL_ERROR "SCR can not be used with MPIIO" )
  endif()
				message (STATUS " MPIIO is defined")
				add_definitions(-DMPIIO)
endif()

##### Optional CP datatypes #####

if(MKL_CP)
  add_definitions(-DMKL_CP)
  message (STATUS "MKL_CP is defined")
  find_path(MKL_INCLUDE_DIR mkl.h $ENV{MKLROOT}/include)
  set(MKL_LIBS $ENV{MKL_LIB})
  include_directories (${MKL_INCLUDE_DIR} ) 
  set(LIBS ${LIBS} ${MKL_LIBS} dl)

  set (MKL_INC_DIR 	${PROJECT_INC_DIR}/cpTypes/cpMkl/)
	set (MKL_INC 	${MKL_INC_DIR}/cpMkl.hpp)
	set (MKL_INC 	${MKL_INC_DIR}/cpMklArray.hpp)
	set (HEADERS 	${HEADERS} ${MKL_INC})
	include_directories (${MKL_INC_DIR})
  set (CPMKL_SRC         ${MKL_INC_DIR}/cpMkl.cpp)
  set (SOURCE ${SOURCE} ${CPMKL_SRC})
  add_definitions(-DMKL_CP)
endif()

if(GHOST_CP)
  message (STATUS "GHOST_CP is defined")
	find_package(GHOST REQUIRED)
	add_definitions(-DGHOST_CP)
	if(GHOST_FOUND)
		message (STATUS "GHOST pacage is found ${GHOST_INCLUDE_DIRS} ")
		include_directories(${GHOST_INCLUDE_DIRS})
		set(LIBS ${LIBS} ${GHOST_LIBRARIES})
		message (STATUS "GHOST LIBS ${GHOST_LIBRARIES} ")
	else()
		message (STATUS "GHOST package is NOT found")
	endif()
	
	##### CPGHOST #####
	set (CPGHOST_INC_DIR 	${PROJECT_INC_DIR}/cpTypes/cpGhost/)
	set (CPGHOST_INC 	${CPGHOST_INC_DIR}/cpGhost.hpp)
	set (HEADERS 	${HEADERS} ${CPGHOST_INC})
	include_directories (${CPGHOST_INC_DIR})
  set (CPGHOST_SRC         ${CPGHOST_INC_DIR}/cpGhost.cpp)
  set (SOURCE ${SOURCE} ${CPGHOST_SRC})
  
endif()

if(PHIST_CP)
  find_path(PHIST REQUIRED)
	add_definitions(-DPHIST_CP)
	include_directories ("${PHIST}/include/")
	link_directories ("${PHIST}/lib/")
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)

	##### MVEC #####
	set (CPPHISTMVEC_INC_DIR 	${PROJECT_INC_DIR}/cpTypes/cpPhistMvec/)
	include_directories ("${CPPHISTMVEC_INC_DIR}")
	set (CPPHISTMVEC_INC 	${CPPHISTMVEC_INC_DIR}/cpPhistMvec.h
												${CPPHISTMVEC_INC_DIR}/cpPhistMvec_decl.h
												${CPPHISTMVEC_INC_DIR}/cpPhistMvec_def.hpp)
	set (CPPHISTMVEC_SRC 	${CPPHISTMVEC_INC_DIR}/cpPhistMvec.cpp)

	set ( SOURCE 	${SOURCE} ${CPPHISTMVEC_SRC})
	set ( HEADERS 	${HEADERS} ${CPPHISTMVEC_INC})

	##### SDMAT #####
	set (CPPHISTSDMAT_INC_DIR 	${PROJECT_INC_DIR}/cpTypes/cpPhistSdMat/)
	include_directories ("${CPPHISTSDMAT_INC_DIR}")
	set (CPPHISTSDMAT_INC 	${CPPHISTSDMAT_INC_DIR}/cpPhistSdMat.h
													${CPPHISTSDMAT_INC_DIR}/cpPhistSdMat_decl.h
													${CPPHISTSDMAT_INC_DIR}/cpPhistSdMat_def.hpp)
	set (CPPHISTSDMAT_SRC 	${CPPHISTSDMAT_INC_DIR}/cpPhistSdMat.cpp)

	set ( SOURCE 	${SOURCE} ${CPPHISTSDMAT_SRC})
	set ( HEADERS 	${HEADERS} ${CPPHISTSDMAT_INC})

endif()



add_library(craft SHARED ${SOURCE})
set(LIBRARY_INSTALL_DIR lib)

configure_file("${PROJECT_DIR}/craft-config.cmake.in" "${PROJECT_DIR}/craft-config.cmake" @ONLY)


##### INSTALL TARGET #####
install ( TARGETS craft EXPORT craftLibrary LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR} )
install ( FILES ${PROJECT_DIR}/craft-config.cmake DESTINATION ${LIBRARY_INSTALL_DIR} )
	
set(INCLUDE_INSTALL_DIR include)
set(TOOLS_INSTALL_DIR tools)
#install ( FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR} )
install ( DIRECTORY ${PROJECT_INC_DIR} DESTINATION ${INCLUDE_INSTALL_DIR} )
install ( DIRECTORY ${PROJECT_TOOLS_DIR} DESTINATION ${TOOLS_INSTALL_DIR} USE_SOURCE_PERMISSIONS )

##### UNINSTALL TARGET #####
configure_file(
				    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
						 "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
	    				IMMEDIATE @ONLY)

add_custom_target(uninstall
		COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

