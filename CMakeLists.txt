cmake_minimum_required(VERSION 2.8)

# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, 
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used), Debug, Release, RelWithDebInfo")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used), Debug, Release, RelWithDebInfo")
ENDIF()

message(STATUS "build type: ${CMAKE_BUILD_TYPE}")

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Trilinos")

# Define some custom variables
project(PHIST)
set(PHIST_VERSION_MAJOR 0)
set(PHIST_VERSION_MINOR 5)
set(PHIST_VERSION_PATCH 0)
set(PHIST_VERSION_PATCH ${PHIST_VERSION_PATCH}-dev)
set(PRODUCT_NAME phist-${PHIST_VERSION_MAJOR}.${PHIST_VERSION_MINOR}.${PHIST_VERSION_PATCH})

option(PHIST_ENABLE_MPI "Enable MPI within PHIST" ON)

if (DEFINED LAPACK_LIBS)
  set(LAPACK_LIBS ${LAPACK_LIBS} CACHE STRING "how to link with lapack (and blas)")
else()
  set(LAPACK_LIBS "lapack;blas" CACHE STRING "how to link with lapack (and blas)")
endif()

enable_language (C)
enable_language (CXX)
enable_language (Fortran)

##############################################################################
# compiler-dependent flags 
##############################################################################

set(OUTLEV_DEBUG 4)
set(OUTLEV_RELEASE 3)

if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  message(STATUS "using GNU compiler flags")
  set (CMAKE_C_FLAGS "-fopenmp")
  set (CMAKE_C_FLAGS_RELEASE "-O2 -fopenmp")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g -fopenmp")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -g -fopenmp -Wall -Wextra -Wno-format -Wno-sign-compare")

  set (CMAKE_CXX_FLAGS "-fopenmp")
  set (CMAKE_CXX_FLAGS_RELEASE "-O2 -fopenmp")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -fopenmp")
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -DTESTING -fopenmp -Wall -Wextra -Wno-format -Wno-sign-compare")
else ()
#TODO - add flags for Intel and others
  message ("compiler ID: " ${CXX_COMPILER_ID})
  message ("No custom compiler flags are known.")
endif ()

# append output level to C and C++ flags
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -DPHIST_OUTLEV=${OUTLEV_RELEASE}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DPHIST_OUTLEV=${OUTLEV_RELEASE}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DPHIST_OUTLEV=${OUTLEV_DEBUG}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DPHIST_OUTLEV=${OUTLEV_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DPHIST_OUTLEV=${OUTLEV_RELEASE}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPHIST_OUTLEV=${OUTLEV_DEBUG}")


##############################################################################
# setup MPI
##############################################################################

if (${PHIST_ENABLE_MPI})
#TODO - what is the best way to do this portably?
#  find_package(MPI REQUIRED)
#  if (${MPI_FOUND})
    add_definitions(-DPHIST_HAVE_MPI)
#    set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
#    set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
#    set (CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
#    set (CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
#    set (CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
#    set (CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
#    include_directories (${MPI_INCLUDE_PATH})
#  else()
#    message(WARNING "could not find MPI")
#  endif()
set(MPIEXEC "mpirun")
set(MPIEXEC_NUMPROC_FLAG "-np")
endif()

##############################################################################
# setup which kernel library is being used.
##############################################################################

if (NOT PHIST_KERNEL_LIB)
  set(PHIST_KERNEL_LIB $ENV{PHIST_KERNEL_LIB})
endif()

if (NOT PHIST_KERNEL_LIB)
  message(FATAL_ERROR "PHIST_KERNEL_LIB not set.")
endif()

message(STATUS "looking for " ${PHIST_KERNEL_LIB} " ...")

set(LOAD_TRILI_PKG,"")
if  (PHIST_KERNEL_LIB STREQUAL "epetra")
  list(APPEND LOAD_TRILI_PKG "EpetraExt")
  list(APPEND LOAD_TRILI_PKG "Belos")
  add_definitions(-DPHIST_KERNEL_LIB_EPETRA)
elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")
  list(APPEND LOAD_TRILI_PKG "Tpetra")
  list(APPEND LOAD_TRILI_PKG "Belos")
  add_definitions(-DPHIST_KERNEL_LIB_TPETRA)
elseif (PHIST_KERNEL_LIB STREQUAL "ghost")
  # we already link with ghost because of the queuing system  
  list(APPEND LOAD_TRILI_PKG "Belos")
  add_definitions(-DPHIST_KERNEL_LIB_GHOST)
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set" )
endif ()


# these specific versions of type-generic drivers have to be created
if  (PHIST_KERNEL_LIB STREQUAL "epetra")
set(TG_Types D)
else ()
set(TG_Types S 
             D 
             C 
             Z)
endif ()
set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp)
set (PHIST_TOOLS_SOURCE ${PROJECT_SOURCE_DIR}/src/tools/phist_tools.c
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_rcp_helpers.cpp)
set (PHIST_CORE_SOURCE ${PROJECT_SOURCE_DIR}/src/core/phist_operator.c
                       ${PROJECT_SOURCE_DIR}/src/core/phist_orthog.cpp)
#                       ${PROJECT_SOURCE_DIR}/src/core/phist_taskbuf.c)
set (PHIST_KRYLOV_SOURCE ${PROJECT_SOURCE_DIR}/src/krylov/phist_lanczos.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_arnoldi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_gmres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_bgmres.cpp)
set (PHIST_JADA_SOURCE ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOp.c
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jdqr.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_schur_decomp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_blockjada.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejada.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_mvec_times_sdMat_inplace.cpp)


include_directories(${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB})
include_directories(${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/krylov)
include_directories(${PROJECT_SOURCE_DIR}/src/jada)
include_directories(${PROJECT_SOURCE_DIR}/drivers)

add_library(phist_tools STATIC ${PHIST_TOOLS_SOURCE})
add_library(phist_core STATIC ${PHIST_CORE_SOURCE})
add_library(phist_kernels_${PHIST_KERNEL_LIB} STATIC ${PHIST_KERNEL_SOURCE})
add_library(phist_solvers STATIC ${PHIST_KRYLOV_SOURCE} ${PHIST_JADA_SOURCE})

list(APPEND PHIST_LIBS phist_solvers)
list(APPEND PHIST_LIBS phist_core)
list(APPEND PHIST_LIBS phist_kernels_${PHIST_KERNEL_LIB})
list(APPEND PHIST_LIBS phist_tools)

# ghost is a requirement of phist now because of the task queue
#find_package(Ghost PATHS $ENV{GHOST_HOME})
#message(STATUS "GHOST_Found=${GHOST_Found}")
find_package(GHOST REQUIRED)

if (GHOST_FOUND)
  add_definitions(-DPHIST_HAVE_GHOST)
  include_directories(${GHOST_INCLUDE_DIRS})
  link_directories(${GHOST_LIBRARY_DIRS})
  list(APPEND PHIST_LIBS ${GHOST_LIBRARIES})
else()
  message(WARNING "Ghost is required for the scheduling in phist, without it things may "
                " not be functional.                                                    "
                "Append the ghost installation directory to the environment variable    "
                "PKG_CONFIG_PATH to help cmake find it.                                 ")
endif()

find_package(Likwid PATHS "$ENV{LIKWID_HOME}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" Optional)
if (${LIKWID_Found})
  add_definitions(-DPHIST_HAVE_LIKWID)
  include_directories(${LIKWID_INCLUDE_DIRS})
  link_directories(${LIKWID_LIBRARY_DIRS})
  list(APPEND PHIST_LIBS ${LIKWID_LIBRARIES} -pthread)
  if (${LIKWID_PERFMON})
    message(STATUS "activated Likwid performance monitoring")
    add_definitions(-DLIKWID_PERFMON)
  endif ()
else()
  message(STATUS "continuing without Likwid, you can set the environment variable LIKWID_HOME to help CMake find it")
endif()

if( ${PHIST_TIMEMONITOR} )
  message(STATUS "activated TimeMonitor timing monitoring")
  add_definitions(-DPHIST_TIMEMONITOR)
endif()


foreach(PKG ${LOAD_TRILI_PKG})
  find_package(${PKG} REQUIRED)
  if (${${PKG}_FOUND})
    message(STATUS "Found.")
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    list(APPEND PHIST_LIBS ${${PKG}_LIBRARIES})
    list(APPEND PHIST_LIBS ${${PKG}_TPL_LIBRARIES})
  else()
    message(FATAL_ERROR ${PKG} " was not found")
  endif()

endforeach()

list(APPEND PHIST_LIBS ${LAPACK_LIBS})
list(APPEND PHIST_LIBS ${MPI_CXX_LIBRARIES})

#message(STATUS "all libraries: "${PHIST_LIBS})

set(GCOVR_EXCLUDES -e src/kernels/kernels_noimpl.c
                   -e src/tools/tools.c)

# Custom properties
set(INTEGRATION_BUILD 0 CACHE STRING "Indicates an integration build")
set(PYTHON_INTERPRETER python CACHE STRING "String Path to the Python interpreter")
set(GCOVR_SCRIPT /usr/bin/gcovr CACHE STRING "Path to the gcovr script")
set(CPPCHECK_EXECUTABLE /usr/bin/cppcheck CACHE STRING "Path to the cppcheck executable")
if (${INTEGRATION_BUILD})
    set(GCOVR_OPTIONS -r src ${GCOVR_EXCLUDES} --xml -o coverage.xml CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS src --enable=all --xml 2> cppcheck.xml CACHE STRING "cppcheck options")
else()
    set(GCOVR_OPTIONS -r src ${GCOVR_EXCLUDES} CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS src --enable=all CACHE STRING "cppcheck options")
endif()

# Define includes
include_directories(lib/googletest-1.6 src/ src/tools/ src/kernels/)

# Define sources
set(SOURCES 
  ${PHIST_KERNEL_SOURCE}
  ${PHIST_TOOLS_SOURCE}
  ${PHIST_CORE_SOURCE}
  ${PHIST_KRYLOV_SOURCE}
  ${PHIST_JADA_SOURCE}
)

set(TEST_SOURCES
    lib/googletest-1.6/gmock-gtest-all.cc
    test/mainTest.cpp
    test/tools/MatrixIO.cpp)
    
file(GLOB TestMatricesMM ${CMAKE_SOURCE_DIR}/test/matrices/*.mm)
file(GLOB TestMatricesBIN ${CMAKE_SOURCE_DIR}/test/matrices/*.bin)

foreach (mat ${TestMatricesMM})
  get_filename_component(${mat} mat_name NAME)
  configure_file(${mat} ${CMAKE_CURRENT_BINARY_DIR}/${mat_name} COPYONLY)
  list(APPEND TEST_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${mat_name})
endforeach()
foreach (mat ${TestMatricesBIN})
  get_filename_component(${mat} mat_name NAME)
  configure_file(${mat} ${CMAKE_CURRENT_BINARY_DIR}/${mat_name} COPYONLY)
  list(APPEND TEST_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${mat_name})
endforeach()

#######################
# Define executables  #
#######################

# type generic driver routines are compiled from .c or .cpp files
# which we get from .in files by preprocessing
file(GLOB TG_C_Drivers ${CMAKE_SOURCE_DIR}/drivers/*.c.in )
file(GLOB TG_CXX_Drivers ${CMAKE_SOURCE_DIR}/drivers/*.cpp.in )


foreach (infile ${TG_C_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.c)
    add_executable(${tc}${filebase} ${tc}${filebase}.c)
    foreach(lib ${PHIST_LIBS})
      target_link_libraries(${tc}${filebase} ${lib})
    endforeach()
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(${tc}${filebase} ${tc}${filebase}.cpp)
    foreach(lib ${PHIST_LIBS})
      target_link_libraries(${tc}${filebase} ${lib})
    endforeach()
  endforeach()
endforeach()

# test drivers
if( "${CMAKE_BUILD_TYPE}" MATCHES "Release" )
  set(PHIST_TESTS kernels
                  core
                  jada
                  krylov
                  bench)
  #                sched
else()
  set(PHIST_TESTS kernels
                  core
                  jada
                  krylov)
endif()

foreach(label ${PHIST_TESTS})
  string(TOUPPER ${label} label_uc)

  file(GLOB ${label_uc}_TEST_SOURCES ${CMAKE_SOURCE_DIR}/test/${label}/*.cpp)
  list(APPEND ${label_uc}_TEST_SOURCES ${TEST_SOURCES})

  add_executable(${PRODUCT_NAME}-${label}-test ${${label_uc}_TEST_SOURCES} ${TEST_DEPENDS})
  
  foreach(lib ${PHIST_LIBS})
    target_link_libraries(${PRODUCT_NAME}-${label}-test ${lib})
  endforeach()
endforeach()

# some example drivers which haven't made their way into the drivers/ directory
file(GLOB C_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.c )

foreach (infile ${C_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  add_executable(${filebase} ${infile})
  foreach(lib ${PHIST_LIBS})
    target_link_libraries(${filebase} ${lib})
  endforeach()
endforeach()

# type generic examples routines are compiled from .c or .cpp files
# which we get from .in files by preprocessing
file(GLOB TG_C_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.c.in )
file(GLOB TG_CXX_Examples ${CMAKE_SOURCE_DIR}/examples/*/*.cpp.in )

foreach (infile ${TG_C_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.c)
    add_executable(${tc}${filebase} ${tc}${filebase}.c)
    foreach(lib ${PHIST_LIBS})
      target_link_libraries(${tc}${filebase} ${lib})
    endforeach()
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Examples})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(${tc}${filebase} ${tc}${filebase}.cpp)
    foreach(lib ${PHIST_LIBS})
      target_link_libraries(${tc}${filebase} ${lib})
    endforeach()
  endforeach()
endforeach()

#

# Activate testing
enable_testing()

if (${PHIST_ENABLE_MPI})
  set(PHIST_MPI_TEST_NPROCS 1 2 4 7 8)
else()
  set(PHIST_MPI_TEST_NPROCS "")
endif()

foreach(np ${PHIST_MPI_TEST_NPROCS})
  foreach(label ${PHIST_TESTS})
    set(TEST_COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${PRODUCT_NAME}-${label}-test --gtest_output=xml:test_output_np${np}/${label}/xunit.xml)
    add_test(NAME test-${label}_np${np} COMMAND ${TEST_COMMAND})
  endforeach()
endforeach()

# Adding concurrency tests with valgrind-helgrind.
# This was done by Timo for his bachelor thesis, currently
# the whole tasking thing is not working and the use/portability
# of automated valgrind tests is questionable, so we skip this test for now.
#add_test(concurrency-test: ${CMAKE_SOURCE_DIR}/test/multiple_helgrind_run.sh ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR})

# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure)

# Activate gcov code generation
if (CMAKE_COMPILER_IS_GNUCXX AND "${CMAKE_BUILD_TYPE}" MATCHES "Debug")
    message (STATUS "Coverage support enabled for targets: ${ARG_TARGETS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage -lgcov")
endif()

# Custom targets: doc, cppcheck, coverage, audit
# Adding optinal target "doc" to generate Doxygen documentation
find_package(Doxygen)
if (DOXYGEN_FOUND)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/phist.doxyfile
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# Add optional target "coverage" to genereate code coverage report (console / XML)
if (GCOVR_SCRIPT)
    add_custom_target(coverage
    ${PYTHON_INTERPRETER} ${GCOVR_SCRIPT} ${GCOVR_OPTIONS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# Adding optional "cppcheck" target
if (CPPCHECK_EXECUTABLE)
    add_custom_target(cppcheck
    ${CPPCHECK_EXECUTABLE} ${CPPCHECK_OPTIONS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()

# Audit target: Combines check, cppcheck, coverage
add_custom_target(audit DEPENDS coverage cppcheck)

# Define the distribution package
set(CPACK_PACKAGE_VERSION "${PHIST_VERSION_MAJOR}.${PHIST_VERSION_MINOR}.${PHIST_VERSION_PATCH}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_GENERATOR ZIP)
#install(TARGETS ${PRODUCT_NAME} RUNTIME DESTINATION bin)
install(FILES "${CMAKE_SOURCE_DIR}/AUTHORS.txt" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/CHANGES.txt" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/LICENSE.txt" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/README.txt" DESTINATION .)
install(DIRECTORY "${CMAKE_SOURCE_DIR}/doc/html" DESTINATION doc)
include(CPack)
