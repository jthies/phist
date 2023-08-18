set(REQUIRE_TRILI_PKG,"")
set(OPTIONAL_TRILI_PKG,"")

if (PHIST_USE_SOLVER_TPLS)
  list(APPEND OPTIONAL_TRILI_PKG "Teuchos")
  list(APPEND OPTIONAL_TRILI_PKG "Anasazi")
  list(APPEND OPTIONAL_TRILI_PKG "Belos")
endif()

if  (PHIST_KERNEL_LIB STREQUAL "epetra")
  list(APPEND REQUIRE_TRILI_PKG "Epetra")
  list(APPEND REQUIRE_TRILI_PKG "EpetraExt")
  list(APPEND REQUIRE_TRILI_PKG "Teuchos")

  list(APPEND OPTIONAL_TRILI_PKG "Kokkos")
  if (PHIST_USE_GRAPH_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Isorropia")
    list(APPEND OPTIONAL_TRILI_PKG "Zoltan")
  endif()
  if (PHIST_USE_PRECON_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Ifpack")
    list(APPEND OPTIONAL_TRILI_PKG "ML")
    list(APPEND OPTIONAL_TRILI_PKG "MueLU")
  endif()
elseif  (PHIST_KERNEL_LIB STREQUAL "tpetra")
  list(APPEND REQUIRE_TRILI_PKG "Tpetra")
  list(APPEND REQUIRE_TRILI_PKG "Kokkos")
  list(APPEND REQUIRE_TRILI_PKG "Teuchos")
  if (PHIST_USE_PRECON_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Ifpack2")
    list(APPEND OPTIONAL_TRILI_PKG "Amesos2")
    list(APPEND OPTIONAL_TRILI_PKG "MueLU")
  endif()

elseif  (PHIST_KERNEL_LIB STREQUAL "ghost")

  if (PHIST_USE_GRAPH_TPLS)
    list(APPEND OPTIONAL_TRILI_PKG "Zoltan")
  endif()
  list(APPEND OPTIONAL_TRILI_PKG "Kokkos")
  list(APPEND OPTIONAL_TRILI_PKG "Teuchos")

endif()


###########################################################################################
# check for required Trilinos libraries, note that you have to add the Trilinos lib/cmake #
# dir to your CMAKE_PREFIX_PATH so that CMake finds them automatically.                   #
###########################################################################################

foreach(PKG ${REQUIRE_TRILI_PKG})
  string(TOUPPER ${PKG} PKG_CAPS)
  if (Trilinos_DIR)
    set(${PKG}_DIR "${Trilinos_DIR}/../${PKG}")
    # Note: originally, we required all trilinos packages to be taken from Trilinos_DIR
    # here, but since Trilinos 14.4.0, Kokkos can be compiled outside of Trilinos, and that
    # caused cmake to fail to find it here. So now we first try to find the packages inside
    # Trilinos_DIR, and if they are not there, we look for them anywhere alse.
    find_package(${PKG} PATHS ${${PKG}_DIR} QUIET NO_DEFAULT_PATH)
    if (NOT ${${PKG}_FOUND})
      find_package(${PKG} REQUIRED)
    endif()
  else()
    find_package(${PKG} REQUIRED)
  endif()
  message(STATUS "Look for required Trilinos package ${PKG}...")
  if (${${PKG}_FOUND})
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    set(PHIST_HAVE_${PKG_CAPS} "Enable required Trilinos package ${PKG}" ON)
    mark_as_advanced(${PKG}_DIR)
  else()
    message(FATAL_ERROR "required package ${PKG} was not found")
  endif()

endforeach()

#########################################################################################
# check for and enable optional libs from Trilinos (preconditioners, Belos/Anasazi)     #
#########################################################################################


if (NOT PHIST_USE_TRILINOS_TPLS)
  set(OPTIONAL_TRILI_PKG "")
endif()

foreach(PKG ${OPTIONAL_TRILI_PKG})
  string(TOUPPER ${PKG} PKG_CAPS)
  set(TPL_REPORT "optional Trilinos package ${PKG_CAPS}... ")
  if (Trilinos_DIR)
    set(${PKG}_DIR "${Trilinos_DIR}/../${PKG}")
    find_package(${PKG} PATHS ${${PKG}_DIR} QUIET NO_DEFAULT_PATH QUIET)
  else()
    find_package(${PKG} QUIET)
  endif()

  if (${${PKG}_FOUND})
    list(APPEND TPL_REPORT "found")
    include_directories(${${PKG}_INCLUDE_DIRS})
    include_directories(${${PKG}_TPL_INCLUDE_DIRS})
    link_directories(${${PKG}_LIBRARY_DIRS})
    link_directories(${${PKG}_TPL_LIBRARY_DIRS})
    mark_as_advanced(${PKG}_DIR)
    set(PHIST_HAVE_${PKG_CAPS} "Enable optional Trilinos package ${PKG}" ON)
  else()
    list(APPEND TPL_REPORT "not found")
  endif()
  message(STATUS ${TPL_REPORT})
endforeach()

