# This file contains a list of all third-party libraries
# supported by at least one of the kernel library interfaces.
# It is intended to implement the behavior prescribed by the xSDK
# policy M1, namely to check the variables
# TPL_ENABLE_<name>, TPL_<name>_DIR/LIBRARIES/INCLUDE_DIRS and
# to use the package as specified. We do this by calling find_package
# with corresponding arguments once for ech TPL before the actual calls
# in CMakeLists.txt. For TPLs that provide a <pkg>-config.cmake or <pkg>Config.cmake
# file, we only support setting TPL_<pkg>_DIR because checking wether the given include
# path and library lists are correct is beyond our capabilities. For other packages (e.g. ParMETIS)
# we have our own Find<pkg>.cmake modules and do our best to implement the full M1 beavior, namely to
# honor the coices made via TPL_<pkg>_INCLUDE_DIRS and TPL_<pkg>_LIBRARIES and stop with an error if we
# find they don't work in some way. We certainly recommend against using this mechanism in favor of
# CMAKE_PREFIX_PATH, PKG_CONFIG_PATH, etc., supported by TPL_<pkg>_DIR where needed.
#
set(PHIST_MODULE_TPL_LIST
  ColPack
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

# these packages provide a <pkg>-config.cmake file, and we only support TPL_<pkg>_DIR
foreach (PKG in ${PHIST_CONFIG_TPL_LIST})
  set(TPL_ENABLE_${PKG} ON CACHE BOOL "Try to find and use ${TPL} (if supported by th kernel library).")
  if (TPL_ENABLE_${PKG})
    if (TPL_${PKG}_INCLUDE_DIRS OR TPL_${PKG}_LIBRARIES)
      message(FATAL_ERROR "HandleTPLS: you specified TPL_${PKG}_INCLUDE_DIRS and/or TPL_${PKG}_LIBRARIES, but ${PKG}"
                          "            is a CMake project and should provide a config file. Either add its location "
                          "            to the CMAKE_PREFIX_PATH or use the variable TPL_${PKG}_DIR instead.")
    elseif (TPL_${PKG}_DIR)
      find_package(${PKG} PATHS ${TPL_${PKG}_DIR} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH REQUIRED)
    elseif (TPL_${PKG}_REQUIRED)
      find_package(${PKG} REQUIRED)
    else()
      find_package(${PKG} QUIET)
    endif()
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
      find_package(${PKG} MODULE REQUIRED)
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

