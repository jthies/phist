cmake_minimum_required (VERSION 3.8)
project (SCAMAC C)
set (SCAMAC_VERSION_MAJOR 0)
set (SCAMAC_VERSION_MINOR 8)

## general setup ...
# add the binary tree directory to the search path for include files
include_directories( ${CMAKE_CURRENT_BINARY_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

option(BUILD_SHARED_LIBS "build shared libraries (ON/OFF)" ON)

#Add flags for warning (Clang / GNU)
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_C_COMPILER_ID}" STREQUAL "Clang" OR
    "${CMAKE_C_COMPILER_ID}" STREQUAL "AppleClang")
    #set(C_warnings "-Wall -Wextra -Werror")
    set(C_warnings "-Wall")
    set(CMAKE_C_LINK_FLAGS "-lm")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_warnings}")


# include and library directory
set(SCAMAC_INCDIR include)
set(SCAMAC_LIBDIR lib)

# specify options
OPTION(SCAMAC_USE_64 "Compile library with 64 bit indices" ON)
OPTION(SCAMAC_BUILD_APP "Build ScaMaC toolbox" ON)
OPTION(SCAMAC_USE_PNG "Compile with PNG support" OFF)
OPTION(SCAMAC_USE_BLAS "Compile with BLAS support" OFF)
OPTION(SCAMAC_USE_LAPACK "Compile with LAPACK support" OFF)
OPTION(SCAMAC_USE_EXTERNAL_CBLAS "Compile with external CBLAS" OFF)
OPTION(SCAMAC_USE_EXTERNAL_CLAPACK "Compile with external CLAPACK" OFF)
OPTION(SCAMAC_USE_OPENMP "Compile with OpenMP support" ON)
OPTION(SCAMAC_BUILD_MWE "Build minimal working examples" ON)
#OPTION(SCAMAC_USE_SLEPC "Build SLEPc example(s)" OFF)

# check options for consistency

if (SCAMAC_USE_LAPACK AND NOT SCAMAC_USE_BLAS)
  message(FATAL_ERROR "SCAMAC_USE_LAPACK requires SCAMAC_USE_BLAS")
endif()

# check for types

INCLUDE (CheckTypeSize)
if (SCAMAC_USE_64) 
  set(SCAMAC_INDEX_TYPE_int64 ON)
  CHECK_TYPE_SIZE(int64_t SIZE_INT64)
  if (NOT HAVE_SIZE_INT64)
    message(FATAL_ERROR "int64_t does not exist")
  endif()
else()
  set(SCAMAC_INDEX_TYPE_int32 ON)
  CHECK_TYPE_SIZE(int32_t SIZE_INT32)
  if (NOT HAVE_SIZE_INT32)
    message(FATAL_ERROR "int32_t does not exist")
  endif()
endif()

# check for packages

## put custom Findxxx.cmake modules into ./cmake/Modules
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# put custom Findxxx.cmake modules into .
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")


if (SCAMAC_USE_PNG) 
  find_package(PNG REQUIRED)
  include_directories(${PNG_INCLUDE_DIRS})
endif()

set(SCAMAC_EXTERNAL_CBLAS "")
if (SCAMAC_USE_BLAS)
  find_package(BLAS REQUIRED)
  include_directories(${BLAS_INCLUDE_DIRS})
  if (SCAMAC_USE_EXTERNAL_CBLAS)
    find_package(CBLAS REQUIRED)
    include_directories(${CBLAS_INCLUDE_DIR})
    set(SCAMAC_EXTERNAL_CBLAS "<${CBLAS_INCLUDE_FILE}>")
  endif()
endif()

set(SCAMAC_EXTERNAL_CLAPACK "")
if (SCAMAC_USE_LAPACK)
  find_package(LAPACK REQUIRED)
  include_directories(${LAPACK_INCLUDE_DIRS})
  if (SCAMAC_USE_EXTERNAL_CLAPACK)
    if (CLAPACK_INCLUDE_FILE)
      set(SCAMAC_EXTERNAL_CLAPACK "<${CLAPACK_INCLUDE_FILE}>")
    else()
      set(SCAMAC_EXTERNAL_CLAPACK "<${CBLAS_INCLUDE_FILE}>")
    endif()
  endif()
endif()

if (SCAMAC_USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# leave that to external source ...
#[[
if (SCAMAC_USE_SLEPC)
  if(PETSC_DIR STREQUAL "")
    unset(PETSC_DIR CACHE)
  endif()
  if(PETSC_ARCH STREQUAL "")
    unset(PETSC_ARCH CACHE)
  endif()
  if(SLEPC_DIR STREQUAL "")
    unset(SLEPC_DIR CACHE)
  endif()
  set(PETSC_DIR $ENV{PETSC_DIR} CACHE PATH "PETSc directory")
  set(PETSC_ARCH $ENV{PETSC_ARCH} CACHE PATH "PETSc archive")
  set(SLEPC_DIR $ENV{SLEPC_DIR} CACHE PATH "SLEPc directory")
  set(SCAMAC_INCDIR ${CMAKE_CURRENT_BINARY_DIR}/include)
  set(SCAMAC_LIBDIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
  configure_file(makefile_slepc.in makefile_slepc @ONLY)
  file(COPY ${CMAKE_SOURCE_DIR}/mwe_slepc.c DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(mwe_slepc COMMAND make -f ${CMAKE_CURRENT_BINARY_DIR}/makefile_slepc mwe_slepc)
endif()
]]


# configure headers
configure_file(scamac_config.h.in scamac_config.h)
configure_file(scamac_tools.h.in scamac_tools.h)


# list of files
set(SCAMAC_DOF_FILES scamac_dof_fermions.c scamac_dof_bosons.c scamac_dof_spins.c)
#file(GLOB SCAMAC_DOF_FILES       scamac_dof_*.c)
set(SCAMAC_MATRIX_FILES scamac_matrix_FreeFermionChain.c scamac_matrix_SpinChainXXZ.c scamac_matrix_Tridiagonal.c scamac_matrix_FreeBosonChain.c scamac_matrix_TridiagonalComplex.c scamac_matrix_Anderson.c scamac_matrix_Hubbard.c scamac_matrix_Harmonic.c scamac_matrix_TridiagonalReal.c scamac_matrix_Exciton.c)
set(SCAMAC_MATRIX_HEADERS scamac_matrix_Hubbard.h scamac_matrix_TridiagonalComplex.h scamac_matrix_TridiagonalReal.h scamac_matrix_FreeFermionChain.h scamac_matrix_Tridiagonal.h scamac_matrix_FreeBosonChain.h scamac_matrix_Harmonic.h scamac_matrix_Anderson.h scamac_matrix_SpinChainXXZ.h scamac_matrix_Exciton.h)
set(SCAMAC_WRAPPER_FILES scamac_wrapper_OneFermion.c)
set(SCAMAC_WRAPPER_HEADERS scamac_wrapper_OneFermion.h)
set(SCAMAC_MWE_SERIAL mwe_serial_scale.c mwe_serial_query.c mwe_serial_count.c)
set(SCAMAC_MWE_OPENMP mwe_openmp_statistics.c mwe_openmp_count.c)
set(SCAMAC_MWE_MPI mwe_mpi_count.c)
set(SCAMAC_MWE_SLEPC mwe_slepc.c)





#set(SCAMAC_LIB_FILES scamac_aux.c scamac_safeint.c scamac_rng.c scamac_multidx.c scamac_sparserow.c ${SCAMAC_DOF_FILES} scamac_lut.c ${SCAMAC_MATRIX_FILES} ${SCAMAC_WRAPPER_FILES} scamac_generator.c scamac_sparsemat.c scamac_sparsemat_io.c scamac_statistics.c scamac_collection.c scamac_string.c)
set(SCAMAC_LIB_FILES scamac_aux.c scamac_collection.c scamac_generator.c scamac_lut.c scamac_multidx.c scamac_rng.c scamac_safeint.c scamac_sparserow.c scamac_string.c ${SCAMAC_DOF_FILES} ${SCAMAC_MATRIX_FILES} ${SCAMAC_WRAPPER_FILES})


set(SCAMAC_TOOLKIT_FILES scamac_sparsemat_io.c scamac_sparsemat.c scamac_statistics.c)
if (SCAMAC_USE_BLAS)
  if (NOT SCAMAC_USE_EXTERNAL_CBLAS)
    LIST(APPEND SCAMAC_TOOLKIT_FILES scamac_cblas.c)
  endif()
endif()
if (SCAMAC_USE_LAPACK)
  LIST(APPEND SCAMAC_TOOLKIT_FILES scamac_lanczos.c)
  LIST(APPEND SCAMAC_TOOLKIT_FILES scamac_spectrum.c)
endif()
if (SCAMAC_USE_PNG)
  LIST(APPEND SCAMAC_TOOLKIT_FILES scamac_plot.c)
endif()

# add libraries
if (BUILD_SHARED_LIBS)
  ADD_LIBRARY(scamac SHARED ${SCAMAC_LIB_FILES})
  ADD_LIBRARY(scamactools SHARED ${SCAMAC_TOOLKIT_FILES})
else()
  ADD_LIBRARY(scamac STATIC ${SCAMAC_LIB_FILES})
  ADD_LIBRARY(scamactools STATIC ${SCAMAC_TOOLKIT_FILES})
endif()
if (SCAMAC_BUILD_APP)
  set(SCAMAC_APP_FILES scamac_toolbox_app.c)
  add_executable(scamact ${SCAMAC_APP_FILES})
  target_link_libraries(scamact scamac scamactools m ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${PNG_LIBRARIES})
endif()

if (SCAMAC_BUILD_MWE)
  if (SCAMAC_USE_PNG)
    list(APPEND SCAMAC_MWE_SERIAL mwe_plot_serial.c)
    list(APPEND SCAMAC_MWE_OPENMP mwe_plot_openmp.c)
  endif()
  foreach(file ${SCAMAC_MWE_SERIAL})
    get_filename_component(name_without_extension ${file} NAME_WE)
    add_executable(${name_without_extension} ${file})
    target_link_libraries(${name_without_extension} scamac m)
    list(APPEND SCAMAC_MWE ${name_without_extension})
  endforeach()
  if (SCAMAC_USE_OPENMP) 
    foreach(file ${SCAMAC_MWE_OPENMP})
      get_filename_component(name_without_extension ${file} NAME_WE)
      add_executable(${name_without_extension} ${file})
      target_link_libraries(${name_without_extension} scamac m ${OpenMP_C_FLAGS})
      target_compile_options(${name_without_extension} PRIVATE ${OpenMP_C_FLAGS})
      list(APPEND SCAMAC_MWE ${name_without_extension})
    endforeach()
  endif()
  foreach(file ${SCAMAC_MWE})
  # add scamactools as necessary
    if (${file} MATCHES "statistics")
      target_link_libraries(${file} scamactools)
    endif()
  # add PNG library as necessary
    if (${file} MATCHES "plot")
      target_link_libraries(${file} scamactools ${PNG_LIBRARIES})
    endif()
  endforeach()
  # configure makefiles for external mwe
  configure_file(makefile_mwe_external ${CMAKE_CURRENT_BINARY_DIR}/mwe_external/Makefile       @ONLY)
  configure_file(makefile_slepc.in     ${CMAKE_CURRENT_BINARY_DIR}/mwe_external/makefile_slepc @ONLY)
  configure_file(makefile_mpi.in       ${CMAKE_CURRENT_BINARY_DIR}/mwe_external/makefile_mpi   @ONLY)
  foreach(file ${SCAMAC_MWE_MPI} ${SCAMAC_MWE_SLEPC})
    file(COPY ${file} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/mwe_external)
  endforeach()
  # add_custom_target(mwe_slepc COMMAND make -f ${CMAKE_CURRENT_BINARY_DIR}/makefile_slepc mwe_slepc)
endif()


### install (we install in the build directory)


set_target_properties(scamac PROPERTIES PUBLIC_HEADER "scamac.h;${CMAKE_CURRENT_BINARY_DIR}/scamac_config.h;scamac_include.h;scamac_generator.h;scamac_collection.h")
set_target_properties(scamactools PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_BINARY_DIR}/scamac_tools.h;scamac_statistics.h;scamac_sparsemat.h;scamac_sparsemat_io.h;scamac_plot.h;scamac_lanczos.h;scamac_spectrum.h;scamac_omp.h")

INSTALL(TARGETS scamac scamactools
        ARCHIVE DESTINATION ${SCAMAC_LIBDIR}
        LIBRARY DESTINATION ${SCAMAC_LIBDIR}
        PUBLIC_HEADER DESTINATION ${SCAMAC_INCDIR}
)

INSTALL(TARGETS scamact DESTINATION bin)
#INSTALL(TARGETS ${SCAMAC_MWE} DESTINATION scamac_mwe)

### tests

include(CTest)

# silly tests

add_test (TestRuns scamact)

add_test (TestSimple scamact list) 
set_tests_properties(TestSimple PROPERTIES PASS_REGULAR_EXPRESSION "Hubbard")


add_executable(TestRandom scamac_test_random.c)
target_link_libraries(TestRandom scamac)
add_test(TestRandom TestRandom)


