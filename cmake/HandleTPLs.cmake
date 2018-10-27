# This file contains a list of all third-party libraries
# supported by at least one of the kernel library interfaces.
# It is intended to implement the behavior prescribed by the xSDK
# policy M1, namely to check the variables
# TPL_ENABLE_<name>, TPL_<name>_DIR/LIBRARIES/INCLUDE_DIRS and
# to use the package as specified. We do this by calling find_package
# repeatedly and makinng sure the given paths are the ones CMake finds, too.
#
set(PHIST_MODULE_TPL_LIST
  ColPack
  MPACK
  MAGMA
  ParMETIS
  PETSC
  )

set(PHIST_CONFIG_TPL_LIST
  Eigen3
  Trilinos
  )
set(TPL_Eigen3_REQUIRED ${PHIST_KERNEL_LIB_EIGEN})
set(TPL_MAGMA_REQUIRED ${PHIST_KERNEL_LIB_MAGMA})
set(TPL_PETSC_REQUIRED ${PHIST_KERNEL_LIB_PETSC})
if (PHIST_KERNEL_LIB_EPETRA OR PHIST_KERNEL_LIB_TPETRA)
  set(TPL_Trilinos_REQUIRED ON)
endif()

# these packages provide a <pkg>-config.cmake file, but mostly
# those modules will call find_package internally.
foreach (PKG in ${PHIST_CONFIG_TPL_LIST})
  set(TPL_ENABLE_${PKG} ON CACHE BOOL "Try to find and use ${TPL} (if supported by th kernel library).")
  if (TPL_ENABLE_${PKG})
    if (TPL_${PKG}_INCLUDE_DIRS)
      # this is just to make sure the package is actually installed in that location:
      find_package(${PKG} MODULE HINTS ${TPL_${PKG}_INCLUDE_DIRS} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH REQUIRED)
      # override whatever cmake found with whatever the user specified, and hope for the best...
      set(${PKG}_INCLUDE_DIRS "${TPL_${PKG}_INCLUDE_DIRS}")
      set(${PKG}_LIBRARIES ${TPL_${PKG}_LIBRARIES})
    elseif (TPL_${PKG}_DIR)
      find_package(${PKG} HINTS ${TPL_${PKG}_DIR} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH REQUIRED)
    elseif (TPL_${PKG}_REQUIRED)
      find_package(${PKG} REQUIRED)
    else()
      find_package(${PKG} QUIET)
    endif()
  endif()
endforeach()

# ParMETIS requires METIS, to avoid linking with incompatible libraries, we force the user to set
# both TPL_* options if any
foreach (ITEM in DIR;INCLUDE_DIRS;LIBRARIES)
  if (TPL_ParMETIS_${ITEM})
    if (NOT TPL_METIS_${ITEM})
      message(FATAL_ERROR "If you specify TPL_ParMETIS_${ITEM}, you also have to set TPL_METIS_${ITEM}.")
    endif()
    set(METIS_${ITEM} ${TPL_METIS_${ITEM}})
  endif()
endforeach()

# For these packages we have our own Find<pkg>.cmake module in phist/cmake/.
# We just need to translate the TPL_ options for them to work.
foreach (PKG in ${PHIST_MODULE_TPL_LIST})
  string(TOUPPER ${PKG} PKG_CAPS)
  set(TPL_ENABLE_${PKG} ON CACHE BOOL "Try to find and use ${TPL} (if supported by th kernel library).")
  if (TPL_ENABLE_${PKG})
    if (TPL_${PKG}_INCLUDE_DIRS)
      # override whatever cmake found with whatever the user specified, and hope for the best...
      set(${PKG}_INCLUDE_DIRS "${TPL_${PKG}_INCLUDE_DIRS}")
      set(${PKG_CAPS}_INCLUDE_DIRS "${TPL_${PKG}_INCLUDE_DIRS}")
      set(${PKG}_LIBRARIES ${TPL_${PKG}_LIBRARIES})
      set(${PKG_CAPS}_LIBRARIES ${TPL_${PKG}_LIBRARIES})
      find_package(${PKG} MODULE NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH REQUIRED)
    elseif (TPL_${PKG}_DIR)
      set(${PKG}_DIR ${TPL_${PKG}_DIR})
      set(${PKG_CAPS}_DIR ${TPL_${PKG}_DIR})
      find_package(${PKG} MODULE REQUIRED)
    elseif (TPL_${PKG}_REQUIRED)
      find_package(${PKG} REQUIRED)
    else()
      find_package(${PKG} QUIET)
    endif()
  endif()
endforeach()

# some postprocessing to get the right behavior
if (NOT TPL_ENABLE_Trilinos)
    set(PHIST_USE_TRILINOS_TPLS OFF)
endif()

