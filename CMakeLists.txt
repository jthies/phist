cmake_minimum_required(VERSION 2.8)

set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Trilinos")

# Define some custom variables
project(PHIST)
set(PHIST_VERSION_MAJOR 0)
set(PHIST_VERSION_MINOR 1)
set(PHIST_VERSION_PATCH 0)
set(PHIST_VERSION_PATCH ${PHIST_VERSION_PATCH}-dev)
set(PRODUCT_NAME phist-${PHIST_VERSION_MAJOR}.${PHIST_VERSION_MINOR}.${PHIST_VERSION_PATCH})

option(PHIST_ENABLE_MPI "Enable MPI within PHIST" ON)

enable_language (Fortran)
enable_language (C)
enable_language (CXX)

##############################################################################
# compiler-dependent flags 
##############################################################################

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)


if (C_COMPILER_NAME STREQUAL "gcc")
  # gcc
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
else ()
  message ("CMAKE_C_COMPILER full path: " ${CMAKE_C_COMPILER})
  message ("C compiler: " ${C_COMPILER_NAME})
  message ("No optimized C compiler flags are known, we just try -O2...")
  set (CMAKE_C_FLAGS_RELEASE "-O2")
  set (CMAKE_C_FLAGS_DEBUG   "-O0 -g")
endif ()

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
else ()
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ()

if (CXX_COMPILER_NAME STREQUAL "g++")
  # g++
  set (CMAKE_CXX_FLAGS_RELEASE "-O3")
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")
else ()
  message ("CMAKE_CXX_COMPILER full path: " ${CMAKE_CXX_COMPILER})
  message ("CXX compiler: " ${CXX_COMPILER_NAME})
  message ("No optimized CXX compiler flags are known, we just try -O2...")
  set (CMAKE_CXX_FLAGS_RELEASE "-O2")
  set (CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")
endif ()

##############################################################################
# setup MPI
##############################################################################

if (${PHIST_ENABLE_MPI})
  find_package(MPI REQUIRED)
  if (${MPI_FOUND})
    add_definitions(-DPHIST_HAVE_MPI)
    set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    set (CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
    set (CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
    set (CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
    set (CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
    include_directories (${MPI_INCLUDE_PATH})
  else()
    message(WARNING "could not find MPI")
  endif()
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


set(PHIST_LIBS "")
set(LOAD_TRILI_PKG,"")
if  (PHIST_KERNEL_LIB STREQUAL "epetra")
  list(APPEND LOAD_TRILI_PKG "EpetraExt")
elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")
  list(APPEND LOAD_TRILI_PKG "Tpetra")
elseif (PHIST_KERNEL_LIB STREQUAL "ghost")
  message(FATAL_ERROR "ghost not supported in CMakeLists.txt, yet")
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set" )
endif ()

set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp)
set (PHIST_TOOLS_SOURCE ${PROJECT_SOURCE_DIR}/src/tools/phist_tools.c)

add_library(libphist_kernels_${PHIST_KERNEL_LIB}.a ${PHIST_KERNEL_SOURCE})
add_library(libphist_tools.a ${PHIST_TOOLS_SOURCE})

list(APPEND PHIST_LIBS libphist_kernels_${PHIST_KERNEL_LIB}.a)
list(APPEND PHIST_LIBS libphist_tools.a)

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

list(APPEND PHIST_LIBS ${MPI_CXX_LIBRARIES})

# Custom properties
set(INTEGRATION_BUILD 0 CACHE STRING "Indicates an integration build")
set(PYTHON_INTERPRETER python CACHE STRING "String Path to the Python interpreter")
set(GCOVR_SCRIPT /usr/bin/gcovr CACHE STRING "Path to the gcovr script")
set(CPPCHECK_EXECUTABLE /usr/bin/cppcheck CACHE STRING "Path to the cppcheck executable")
if (${INTEGRATION_BUILD})
    set(GCOVR_OPTIONS -r src --xml -o coverage.xml CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS src --enable=all --xml 2> cppcheck.xml CACHE STRING "cppcheck options")
else()
    set(GCOVR_OPTIONS -r src CACHE STRING "String gcovr options")
    set(CPPCHECK_OPTIONS src --enable=all CACHE STRING "cppcheck options")
endif()

# Define includes
include_directories(lib/googletest-1.6 src/ src/tools/ src/kernels/)

# Define sources
set(SOURCES 
    src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp
    examples/kernels/main_kernel_usage.c
)
set(KERNEL_TEST_SOURCES
    lib/googletest-1.6/gmock-gtest-all.cc
    test/mainTest.cpp
    test/kernels/CommTest.cpp
    test/kernels/MvecTest.cpp
    test/kernels/CrsMatTest.cpp
#    test/kernels/MapTest.cpp
#    test/kernels/SdMatTest.cpp
)
set(KERNEL_TEST_SOURCES ${KERNEL_TEST_SOURCES} ${SOURCES} ${GTEST_SOURCES} ${GMOCK_SOURCES})
list(REMOVE_ITEM KERNEL_TEST_SOURCES "examples/kernels/main_kernel_usage.c")

# Define executables
add_executable(${PRODUCT_NAME}-kernel-test ${KERNEL_TEST_SOURCES})
add_executable(main_kernel_usage ${PROJECT_SOURCE_DIR}/examples/kernels/main_kernel_usage.c)


# Activate testing
enable_testing()
add_test(kernel-tests ${PRODUCT_NAME}-kernel-test --gtest_output=xml:${CMAKE_SOURCE_DIR}/xunit.xml)
# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure)

foreach(lib ${PHIST_LIBS})
#  message("adding lib '" ${lib} "'")
  target_link_libraries(${PRODUCT_NAME}-kernel-test ${lib})
  target_link_libraries(main_kernel_usage ${lib})
endforeach()

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
