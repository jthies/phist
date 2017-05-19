cmake_minimum_required (VERSION 2.6)
project (ESMAC)
set (ESMAC_VERSION_MAJOR 0)
set (ESMAC_VERSION_MINOR 5)


OPTION(ESMAC_USE_PNG "Compile with PNG support" ON)
OPTION(ESMAC_USE_BLAS "Compile with BLAS support" ON)
OPTION(ESMAC_USE_LAPACK "Compile with LAPACK support" ON)
OPTION(ESMAC_USE_64 "Compile with 64 bit indices" ON)

if (ESMAC_USE_LAPACK AND NOT ESMAC_USE_BLAS)
  message(FATAL_ERROR "ESMAC_USE_LAPACK requires ESMAC_USE_BLAS")
endif()

INCLUDE (CheckTypeSize)
if (ESMAC_USE_64) 
  set(ESMAC_INDEX_TYPE_int64 ON)
  CHECK_TYPE_SIZE(int64_t SIZE_INT64)
  if (NOT HAVE_SIZE_INT64)
    message(FATAL_ERROR "int64_t does not exist")
  endif()
else()
  set(ESMAC_INDEX_TYPE_int32 ON)
  CHECK_TYPE_SIZE(int32_t SIZE_INT32)
  if (NOT HAVE_SIZE_INT32)
    message(FATAL_ERROR "int32_t does not exist")
  endif()
endif()

configure_file(esmac_config.h.in esmac_config.h)

# add the binary tree directory to the search path for include files
include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set (CMAKE_C_FLAGS "-Wall")


# put custom Findxxx.cmake modules into ./cmake/Modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


if (ESMAC_USE_PNG OR ESMAC_USE_BLAS OR ESMAC_USE_LAPACK)
  find_package(PkgConfig)
  set( ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:/usr/local/lib/pkgconfig:/opt/X11/lib/pkgconfig")
endif() 

if (ESMAC_USE_PNG) 
  find_package(PNG)
  include_directories(${PNG_INCLUDE_DIRS})
endif()

if (ESMAC_USE_BLAS)
  find_package(GSL)
  include_directories(${GSL_INCLUDE_DIRS})
  find_package(BLAS)
  include_directories(${BLAS_INCLUDE_DIRS})
endif()

if (ESMAC_USE_LAPACK)
  find_package(LAPACK)
  include_directories(${LAPACK_INCLUDE_DIRS})
endif()

# be warned about the risks of globbing
file (GLOB ESMAC_MATRIX_FILES esmac_matrix_*.c)
file(GLOB ESMAC_MATRIX_HEADERS esmac_matrix_*.h)


set(ESMAC_DOF_FILES esmac_dof_fermions.c esmac_dof_bosons.c esmac_dof_spins.c)

set(ESMAC_LIB_FILES esmac_aux.c esmac_rng.c esmac_multidx.c esmac_sparserow.c esmac_qstate.c esmac_mstate.c ${ESMAC_DOF_FILES} esmac_lut.c ${ESMAC_MATRIX_FILES} esmac_generator.c esmac_sparsemat.c esmac_sparsemat_io.c esmac_statistics.c esmac_collection.c)

if (ESMAC_USE_PNG)
  LIST(APPEND ESMAC_LIB_FILES esmac_plot.c)
endif()

if (ESMAC_USE_BLAS)
endif()

if (ESMAC_USE_LAPACK)
  LIST(APPEND ESMAC_LIB_FILES esmac_lanczos.c)
  LIST(APPEND ESMAC_LIB_FILES esmac_spectrum.c)
endif()

# add our library
ADD_LIBRARY(esmac STATIC ${ESMAC_LIB_FILES})

add_executable(esmacapp esmacapp.c)
target_link_libraries(esmacapp m ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${PNG_LIBRARIES} esmac)

find_package(OpenMP)
if (NOT OPENMP_FOUND)
  message(WARNING "OpenMP not found -- will not build minimal working example (mwe)")
else()
  add_executable(mwe mwe.c)
  target_link_libraries(mwe m esmac ${OpenMP_C_FLAGS})
  target_compile_options(mwe PRIVATE ${OpenMP_C_FLAGS})
endif()


### install

INSTALL(TARGETS esmac DESTINATION ${CMAKE_BINARY_DIR}/lib)
INSTALL(TARGETS esmacapp DESTINATION ${CMAKE_BINARY_DIR}/app)

### tests

include(CTest)

add_test (TestRuns esmacapp)

add_test (TestSimple esmacapp list) 
set_tests_properties(TestSimple PROPERTIES PASS_REGULAR_EXPRESSION "Hubbard")


