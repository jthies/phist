# This file contains lists of all third-party libraries
# supported by at least one of the kernel library interfaces.
# It is intended to implement the behavior prescribed by the xSDK
# policy M1, namely to check the variables
# TPL_ENABLE_<name>, TPL_<name>_DIR/LIBRARIES/INCLUDE_DIRS and
# to use the package as specified. We do this by calling find_package
# with corresponding arguments once for ech TPL before the actual calls
# in CMakeLists.txt. We implement three strategies depending on the TPL at hand:
#
# 1. For TPLs that provide a <pkg>-config.cmake or <pkg>Config.cmake
# file, we only support setting TPL_<pkg>_DIR because checking wether the given include
# path and library lists are correct is beyond our capabilities.
# 2. For packages that provide a <pkg>.pc (pkg-config) file, one may set TPL_<pkg>_DIR to the <prefix>/lib/pkgconfig/
# location but do not accept TPL_<pkg>_INCLUDE_DIRS or TPL_<pkg>_LIBRARIES.
# 3. For other packages (e.g. ParMETIS and ColPack)
# we have our own Find<pkg>.cmake modules and do our best to implement the full M1 beavior, namely to
# honor the coices made via TPL_<pkg>_INCLUDE_DIRS and TPL_<pkg>_LIBRARIES and stop with an error if we
# find they don't work in some way.
#
# For categories 1 and 2 we recommennd using CMAKE_PREFIX_PATH and/or PKG_CONFIG_PATH instead of using the TPL_* mechanism.
#
set(PHIST_MODULE_TPL_LIST
  ColPack
  ParMETIS
  )

set(PHIST_PKG_CONFIG_TPL_LIST
  MAGMA
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

# 1. these packages provide a <pkg>-config.cmake file, and we only support TPL_<pkg>_DIR
foreach (PKG in ${PHIST_CONFIG_TPL_LIST})
  set(TPL_ENABLE_${PKG} ON CACHE BOOL "Try to find and use ${TPL} (if supported by the kernel library).")
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

# 2. packages that can be found using pkg-config (a .pc file). We have Find*.cmake
# files for these, but essentially they just consist of two lines. We only support the TPL_<pkg>_DIR for these.
foreach (PKG in ${PHIST_PKG_CONFIG_TPL_LIST})
  set(TPL_ENABLE_${PKG} ON CACHE BOOL "Try to find and use ${TPL} (if supported by the kernel library).")
  if (TPL_ENABLE_${PKG})
    if (TPL_${PKG}_INCLUDE_DIRS OR TPL_${PKG}_LIBRARIES)
      message(FATAL_ERROR "HandleTPLS: you specified TPL_${PKG}_INCLUDE_DIRS and/or TPL_${PKG}_LIBRARIES, but ${PKG}"
                          "            should provide a pkg-config file (<pkg>.pc). Either add its location "
                          "            to the CMAKE_PREFIX_PATH or PKG_CONFIG_PATH, or set the variable TPL_${PKG}_DIR "
                          "            to the path of the .pc file.")
    elseif (TPL_${PKG}_DIR)
      set(PKG_CONFIG_PATH_OLD $ENV{PKG_CONFIG_PATH})
      set(ENV{PKG_CONFIG_PATH} ${TPL_${PKG}_DIR})
      find_package(${PKG} MODULE REQUIRED)
      set(ENV{PKG_CONFIG_PATH} ${PKG_CONFIG_PATH_OLD})
    elseif (TPL_${PKG}_REQUIRED)
      find_package(${PKG} REQUIRED)
    else()
      find_package(${PKG} QUIET)
    endif()
  endif()
endforeach()

# 3. For these packages we have our own Find<pkg>.cmake module in phist/cmake/, respecting the TPL_* options
foreach (PKG in ${PHIST_MODULE_TPL_LIST})
  string(TOUPPER ${PKG} PKG_CAPS)
  set(TPL_ENABLE_${PKG} ON CACHE BOOL "Try to find and use ${TPL} (if supported by the kernel library).")
  if (TPL_ENABLE_${PKG})
    if (TPL_${PKG}_INCLUDE_DIRS)
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

