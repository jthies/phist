cmake_minimum_required(VERSION 2.8.10)

list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Trilinos")

# Define some custom variables
project(PHIST)
set(PHIST_VERSION_MAJOR 0)
set(PHIST_VERSION_MINOR 7)
set(PHIST_VERSION_PATCH 0)
#set(PHIST_VERSION_PATCH ${PHIST_VERSION_PATCH}-beta)
set(PRODUCT_NAME phist-${PHIST_VERSION_MAJOR}.${PHIST_VERSION_MINOR}.${PHIST_VERSION_PATCH})

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug, Release" FORCE)
else()
  set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: Debug, Release" FORCE)
endif()
option(PHIST_ENABLE_MPI "Enable MPI within PHIST" ON)
option(PHIST_ENABLE_SP "Enable single precision functions and drivers" OFF)
if (DEFINED LAPACK_LIBS)
  set(LAPACK_LIBS ${LAPACK_LIBS} CACHE STRING "how to link with lapack (and blas)")
elseif(NOT CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(LAPACK_LIBS "lapack;blas" CACHE STRING "how to link with lapack (and blas)")
endif()

enable_language (C)
enable_language (CXX)
enable_language (Fortran)

# output behavior
set(PHIST_OUTLEV 3 CACHE STRING "PHIST verbosity (default 3, 'VERBOSE')")

##############################################################################
# compiler-dependent flags 
##############################################################################


if (CMAKE_C_COMPILER_ID STREQUAL "GNU")

  message(STATUS "using GNU compiler flags")
  set (CMAKE_C_FLAGS "-fopenmp -std=c99")
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3 -ffast-math -march=native")
  set (CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS} -O2 -g")
  if (CMAKE_C_COMPILER_VERSION VERSION_LESS 4.8)
    set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS} -dH -O0 -g -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter")
  else ()
    set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS} -dH -O0 -g -fstack-protector-all -DTESTING -fopenmp -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-local-typedefs")
  endif ()

  set (CMAKE_CXX_FLAGS "-fopenmp -std=c++0x")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2")
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O2")
  if (CMAKE_C_COMPILER_VERSION VERSION_LESS 4.8)
    set (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -dH -O0 -g -DTESTING -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter")
  else ()
    set (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -dH -O0 -g -fstack-protector-all -DTESTING -fopenmp -Wall -Wextra -Wno-format -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-unused-local-typedefs")
  endif ()

  set (CMAKE_Fortran_FLAGS "-fopenmp -cpp -ffree-line-length-none")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -march=native -O3 -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -dH -O0 -g -DTESTING -Wall -fcheck=all -fstack-protector-all")

elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")

  message(STATUS "using INTEL compiler flags")
  set (CMAKE_C_FLAGS "-openmp -Wno-unused-variable -mkl -std=c99")
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} -O3 -xHOST -no-prec-div")
  set (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS} -check=uninit,stack -fstack-protector -w2 -debug -traceback -DTESTING")

  set (CMAKE_CXX_FLAGS "-openmp -Wno-unused-variable -mkl -std=c++0x")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -xHOST -no-prec-div")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -check=uninit,stack -fstack-protector -w2 -debug -traceback -DTESTING")

  set (CMAKE_Fortran_FLAGS "-openmp -fpp -mkl")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -xHOST -no-prec-div")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -check all -warn all -debug -traceback -DTESTING")

  # for ipo
  # melven: probably not useful or even dangerous?
  #  kernel routines should probably best optimized individually;
  #  for fortran: I'm passing arrays as scalars in fortran (call by reference),
  #               because the kernel routines would be too complex otherwise
  #               (and they should be fast and not complex!).
  #               The compiler best shouldn't be able to know about this...
  #set_property(DIRECTORY . PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

else()

  message ("compiler ID: " ${CMAKE_C_COMPILER_ID})
  message ("No custom compiler flags are known.")

endif ()

# try to detect if we can use C++11 Lambdas
include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
  int main(int argc, char** argv)
  {
    auto lfunc = [&](){ int _argc = argc; };
    return 0;
  }
" PHIST_HAVE_CXX11_LAMBDAS
)
##############################################################################
# set enable precision?
##############################################################################
if (${PHIST_ENABLE_SP})
  set(PHIST_HAVE_SP 1)
endif ()
##############################################################################
# setup MPI
##############################################################################

if (${PHIST_ENABLE_MPI})
  set(PHIST_HAVE_MPI 1)
#TODO - what is the best way to do this portably?
  find_package(MPI REQUIRED)
  if (${MPI_FOUND})
    set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
    set (CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})
    set (CMAKE_C_LINK_FLAGS ${CMAKE_C_LINK_FLAGS} ${MPI_C_LINK_FLAGS})
    set (CMAKE_Fortran_COMPILE_FLAGS ${CMAKE_Fortran_COMPILE_FLAGS} ${MPI_Fortran_COMPILE_FLAGS})
    set (CMAKE_Fortran_LINK_FLAGS ${CMAKE_Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS})
    include_directories (${MPI_INCLUDE_PATH})
  else()
    message(WARNING "could not find MPI. Trying to compile anyway, presuming the compiler/linker knows where to find it.")
  endif()
set(MPIEXEC "mpirun" CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")
endif()

##############################################################################
# setup which kernel library is being used.
##############################################################################

if("${PHIST_KERNEL_LIB}" STREQUAL "")
  set(PHIST_KERNEL_LIB $ENV{PHIST_KERNEL_LIB} CACHE STRING "Select underlying kernel libraries for phist (epetra|tpetra|ghost|fortran)")
endif()

message(STATUS "looking for " ${PHIST_KERNEL_LIB} " ...")

set(LOAD_TRILI_PKG,"")
if  (PHIST_KERNEL_LIB STREQUAL "epetra")
  list(APPEND LOAD_TRILI_PKG "EpetraExt")
  set(PHIST_KERNEL_LIB_EPETRA 1)
elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")
  list(APPEND LOAD_TRILI_PKG "Tpetra")
  set(PHIST_KERNEL_LIB_TPETRA 1)
elseif (PHIST_KERNEL_LIB STREQUAL "ghost")
  # we already link with ghost because of the queuing system  
  set(PHIST_KERNEL_LIB_GHOST 1)
  set(PHIST_MVECS_ROW_MAJOR 1 CACHE STRING "wether multi-vectors should be stored in row-major ordering")
  set(PHIST_SDMATS_ROW_MAJOR 0 CACHE STRING "store sdMats in row-major order (not supported by our algorithms..." )
elseif (PHIST_KERNEL_LIB STREQUAL "fortran")
  set(PHIST_KERNEL_LIB_FORTRAN 1)
  set(PHIST_MVECS_ROW_MAJOR 1)
else ()
  message( FATAL_ERROR "PHIST_KERNEL_LIB not set or not recognized" )
endif ()

# test for Belos (Trilinos library of iterative solvers).
# We recommend linking with Belos because we use its Tsqr
# interface. If Belos is not found, important functionality
# may not work, in particular the mvec_QR function is currently
# implemented via Belos by all kernel libs except "Fortran"

find_package(Belos)

if (${Belos_FOUND})
  list(APPEND LOAD_TRILI_PKG "Belos")
  set(PHIST_HAVE_BELOS 1)
else()
  message(STATUS "Belos not found, some functionality will be disabled and return -99 ('not implemented')")
endif()


# these specific versions of type-generic drivers have to be created
if  (PHIST_KERNEL_LIB STREQUAL "epetra" OR PHIST_KERNEL_LIB STREQUAL "fortran")
set(TG_Types D)
else ()
  if (PHIST_ENABLE_SP)
    set(TG_Types S 
                 D 
                 C 
                 Z)
  else ()
    set(TG_Types D
                 Z)
  endif ()
endif ()

if (PHIST_KERNEL_LIB STREQUAL "fortran")
  set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/env_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/map_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/mrgrnk.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/crsmat_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/parmetis.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/colpack.cpp
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/mvec_module.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/nrm2_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/axpy_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/axpy_kernels_nt.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/gather_scatter_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/dot_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/gemm_kernels_sB.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/gemm_kernels_sC.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/spmvm_kernels_nt.c
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/spmvm_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/kacz_kernels.f90
                          ${PROJECT_SOURCE_DIR}/src/kernels/fortran/sdmat_module.f90)
else()
  set (PHIST_KERNEL_SOURCE ${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/kernels.cpp)
endif()
set (PHIST_TOOLS_SOURCE ${PROJECT_SOURCE_DIR}/src/tools/phist_tools.c
                        ${PROJECT_SOURCE_DIR}/src/tools/phist_rcp_helpers.cpp)
set (PHIST_CORE_SOURCE ${PROJECT_SOURCE_DIR}/src/core/phist_operator.cpp
                       ${PROJECT_SOURCE_DIR}/src/core/phist_orthog.cpp)
set (PHIST_KRYLOV_SOURCE ${PROJECT_SOURCE_DIR}/src/krylov/phist_simple_arnoldi.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_pgmres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_pminres.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_belos.cpp
                         ${PROJECT_SOURCE_DIR}/src/krylov/phist_carp_cg.cpp)
set (PHIST_JADA_SOURCE   ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaOpts.c
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jdqr.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_schur_decomp.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_subspacejada.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_transform_searchspace.cpp
                         ${PROJECT_SOURCE_DIR}/src/jada/phist_jadaCorrectionSolver.cpp)
set (PHIST_FEAST_SOURCE  ${PROJECT_SOURCE_DIR}/src/feast/phist_feastCorrectionSolver.cpp)


include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB})
include_directories(${PROJECT_SOURCE_DIR}/src/tools)
include_directories(${PROJECT_SOURCE_DIR}/src/core)
include_directories(${PROJECT_SOURCE_DIR}/src/krylov)
include_directories(${PROJECT_SOURCE_DIR}/src/jada)
include_directories(${PROJECT_SOURCE_DIR}/src/feast)
include_directories(${PROJECT_SOURCE_DIR}/drivers)

add_library(phist_tools SHARED ${PHIST_TOOLS_SOURCE})
add_library(phist_core SHARED ${PHIST_CORE_SOURCE})
add_library(phist_kernels_${PHIST_KERNEL_LIB} SHARED ${PHIST_KERNEL_SOURCE})
add_library(phist_solvers SHARED ${PHIST_KRYLOV_SOURCE} 
                                 ${PHIST_JADA_SOURCE}
                                 ${PHIST_FEAST_SOURCE})

list(APPEND PHIST_LIBS phist_solvers)
list(APPEND PHIST_LIBS phist_core)
list(APPEND PHIST_LIBS phist_kernels_${PHIST_KERNEL_LIB})
list(APPEND PHIST_LIBS phist_tools)

if(PHIST_KERNEL_LIB STREQUAL "ghost")
  find_package(GHOST REQUIRED CONFIG)
else()
  find_package(GHOST CONFIG QUIET)
endif()
if (GHOST_FOUND)
  message(STATUS "Found GHOST ${GHOST_VERSION}")
  set(PHIST_HAVE_GHOST 1)
  include_directories(${GHOST_INCLUDE_DIRS})
  link_directories(${GHOST_LIBRARY_DIRS})
  list(APPEND PHIST_LIBS ${GHOST_LIBRARIES})
else()
  message(STATUS "GHOST not found, some functionality will be disabled and return -99 ('not implemented')")
endif()

# load parmetis for reordering fortran
if(PHIST_KERNEL_LIB STREQUAL "fortran")
  
  find_package(ParMETIS)
  if(PARMETIS_FOUND)
    set(PHIST_HAVE_PARMETIS On)
    include_directories(${PARMETIS_INCLUDE_DIRS})
    link_directories(${PARMETIS_LIBRARY_DIRS})
    list(APPEND PHIST_LIBS ${PARMETIS_LIBRARIES})
    message(STATUS "ParMETIS found, enabled matrix reordering!")
  else()
    message(STATUS "ParMETIS not found, no matrix reordering!")
  endif()

#  find_package(ColPack)
  if (COLPACK_DIR)
    set(COLPACK_INCLUDE_DIRS ${COLPACK_DIR}/include)
    set(COLPACK_LIBRARIES ${COLPACK_DIR}/lib/libColPack.so)
    set(COLPACK_LIBRARY_DIRS ${COLPACK_DIR}/lib)
    set(COLPACK_FOUND 1)
  endif()

  if(COLPACK_FOUND)
    set(PHIST_HAVE_COLPACK 1)
    include_directories(${COLPACK_INCLUDE_DIRS})
    link_directories(${COLPACK_LIBRARY_DIRS})
    list(APPEND PHIST_LIBS ${COLPACK_LIBRARIES})
    message(STATUS "ColPack found, enabled matrix coloring!")
  else()
    message(STATUS "ColPack not found, no matrix coloring!")
  endif()

endif()

find_package(Likwid PATHS "$ENV{LIKWID_HOME}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/" Optional)
if (${LIKWID_Found})
  set(PHIST_HAVE_LIKWID 1)
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

set(PHIST_TIMEMONITOR On CACHE BOOL "Gather detailed function-level timings")
if( ${PHIST_TIMEMONITOR} )
  message(STATUS "activated TimeMonitor timing monitoring")
endif()
set(PHIST_TIMEMONITOR_PERLINE Off CACHE BOOL "Gather timing information for each line with a PHIST_CHK_* macro!")

configure_file ("${PROJECT_SOURCE_DIR}/src/phist_config.h.in" "${PROJECT_BINARY_DIR}/phist_config.h")

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
  ${PHIST_FEAST_SOURCE}
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
    add_executable(phist_${tc}${filebase} ${tc}${filebase}.c)
    foreach(lib ${PHIST_LIBS})
      target_link_libraries(phist_${tc}${filebase} ${lib})
    endforeach()
    install(TARGETS phist_${tc}${filebase} RUNTIME DESTINATION bin)
  endforeach()
endforeach()

foreach (infile ${TG_CXX_Drivers})
  get_filename_component(filebase ${infile} NAME_WE)
  foreach (tc ${TG_Types})
    STRING(TOLOWER \"phist_gen_${tc}.h\" PHIST_TG_HFILE)
    configure_file(${infile} ${CMAKE_CURRENT_BINARY_DIR}/${tc}${filebase}.cpp)
    add_executable(phist_${tc}${filebase} ${tc}${filebase}.cpp)
    foreach(lib ${PHIST_LIBS})
      target_link_libraries(phist_${tc}${filebase} ${lib})
    endforeach()
    install(TARGETS phist_${tc}${filebase} RUNTIME DESTINATION bin)
  endforeach()
endforeach()

# example drivers implemented in ghost
if  (PHIST_HAVE_GHOST)
# TODO - make this somehow general, or probably we should
#       actually install phist instead and push the spinChain
#       driver to the essex/examples repo...
  find_package(ESSEX-PHYSICS REQUIRED CONFIG)
  message(STATUS "Found ESSEX-Physics ${ESSEX-PHYSICS_VERSION}")

  include_directories(${ESSEX-PHYSICS_INCLUDE_DIRS})

  add_executable(phist_spinChain ${PROJECT_SOURCE_DIR}/drivers/spinChain.c)
  add_executable(phist_spinChain_subspacejada ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada.cpp ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada_helper.c)
  add_executable(phist_graphene_carp ${PROJECT_SOURCE_DIR}/drivers/graphene_carp.c ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada_helper.c)
#  add_executable(phist_graphene_carp_as_ghost_task ${PROJECT_SOURCE_DIR}/drivers/graphene_carp_as_ghost_task.c ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada_helper.c)
  target_link_libraries(phist_spinChain ${ESSEX-PHYSICS_LIBRARIES})
  target_link_libraries(phist_spinChain_subspacejada ${ESSEX-PHYSICS_LIBRARIES})
  target_link_libraries(phist_graphene_carp ${ESSEX-PHYSICS_LIBRARIES})
#  target_link_libraries(phist_graphene_carp_as_ghost_task ${ESSEX-PHYSICS_LIBRARIES})
  foreach(lib ${PHIST_LIBS})
    target_link_libraries(phist_spinChain ${lib})
    target_link_libraries(phist_spinChain_subspacejada ${lib})
    target_link_libraries(phist_graphene_carp ${lib})
#    target_link_libraries(phist_graphene_carp_as_ghost_task ${lib})
  endforeach()
  install(TARGETS phist_spinChain RUNTIME DESTINATION bin)
  install(TARGETS phist_spinChain_subspacejada RUNTIME DESTINATION bin)
elseif(PHIST_KERNEL_LIB STREQUAL "fortran")
# even without the ghost/physics libs we include certain matfuncs for the
# test cases SpinSZ and Graphene for testing purposes.
  include_directories(${PROJECT_SOURCE_DIR}/drivers/matfuncs)
  add_executable(phist_spinChain              ${PROJECT_SOURCE_DIR}/drivers/spinChain.c
                                        ${PROJECT_SOURCE_DIR}/drivers/matfuncs/SpinChainSZ.c)
  add_executable(phist_spinChain_subspacejada ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada.cpp
                                        ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada_helper.c
                                        ${PROJECT_SOURCE_DIR}/drivers/matfuncs/SpinChainSZ.c)
  add_executable(phist_graphene_carp ${PROJECT_SOURCE_DIR}/drivers/graphene_carp.c
                                        ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada_helper.c
                                        ${PROJECT_SOURCE_DIR}/drivers/matfuncs/graphene.c)
  add_executable(phist_graphene_carp_as_ghost_task ${PROJECT_SOURCE_DIR}/drivers/graphene_carp_as_ghost_task.c
                                        ${PROJECT_SOURCE_DIR}/drivers/spinChain_subspacejada_helper.c
                                        ${PROJECT_SOURCE_DIR}/drivers/matfuncs/graphene.c)
  foreach(lib ${PHIST_LIBS})
    target_link_libraries(phist_spinChain ${lib})
    target_link_libraries(phist_spinChain_subspacejada ${lib})
    target_link_libraries(phist_graphene_carp ${lib})
    target_link_libraries(phist_graphene_carp_as_ghost_task ${lib})
  endforeach()
  install(TARGETS phist_spinChain RUNTIME DESTINATION bin)
  install(TARGETS phist_spinChain_subspacejada RUNTIME DESTINATION bin)
  install(TARGETS phist_graphene_carp RUNTIME DESTINATION bin)
  install(TARGETS phist_graphene_carp_as_ghost_task RUNTIME DESTINATION bin)
endif()

install(TARGETS phist_kernels_${PHIST_KERNEL_LIB} EXPORT phist-targets LIBRARY DESTINATION lib)
install(TARGETS phist_tools EXPORT phist-targets LIBRARY DESTINATION lib)
install(TARGETS phist_core EXPORT phist-targets LIBRARY DESTINATION lib)
install(TARGETS phist_solvers EXPORT phist-targets LIBRARY DESTINATION lib)

install(FILES ${PROJECT_BINARY_DIR}/phist_config.h DESTINATION include/phist)

file(GLOB PHIST_HEADERS ${CMAKE_SOURCE_DIR}/src/*/*.h ${CMAKE_SOURCE_DIR}/src/kernels/${PHIST_KERNEL_LIB}/phist_typedefs.h)
install(FILES ${PHIST_HEADERS} DESTINATION include/phist)

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
  add_executable(phist_${filebase} ${infile})
  foreach(lib ${PHIST_LIBS})
    target_link_libraries(phist_${filebase} ${lib})
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
    set(TEST_COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ./${PRODUCT_NAME}-${label}-test --gtest_output=xml:test_output_np${np}/${label}/xunit.xml)
    add_test(NAME test-${label}_np${np} COMMAND ${TEST_COMMAND})
    set_tests_properties(test-${label}_np${np} PROPERTIES PASS_REGULAR_EXPRESSION "[. \n]*PASSED[. \n]*")
    set_tests_properties(test-${label}_np${np} PROPERTIES FAIL_REGULAR_EXPRESSION "[. \n]*FAILED[. \n]*")
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
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/phist.doxyfile WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
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
#install(FILES "${CMAKE_SOURCE_DIR}/AUTHORS.txt" DESTINATION .)
#install(FILES "${CMAKE_SOURCE_DIR}/CHANGES.txt" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/LICENSE.txt" DESTINATION .)
install(FILES "${CMAKE_SOURCE_DIR}/README.txt" DESTINATION .)
#if (DOXYGEN_FOUND)
#  install(DIRECTORY "${CMAKE_SOURCE_DIR}/doc/html" DESTINATION doc)
#endif()
include(CPack)
