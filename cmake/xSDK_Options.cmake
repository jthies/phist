option(USE_XSDK_DEFAULTS "implement the default behavior suggested for xSDK packages (see https://xsdk.info)" OFF)

option(XSDK_ENABLE_CXX
        "generate C++ bindings for PHIST"
        ON)

option(XSDK_ENABLE_Fortran
        "generate Fortran 2003 bindings for PHIST"
        ON)

if (USE_XSDK_DEFAULTS)
  set (XSDK_INDEX_SIZE 64 CACHE INTEGER "data type for global indices (32/64), default 64.")
  #otherwise, leave it to the kernel library to decide.
endif()

if (XSDK_INDEX_SIZE EQUAL 32)
  set(PHIST_FORCE_32BIT_GIDX_DEFAULT ON)
elseif (XSDK_INDEX_SIZE EQUAL 64)
  set(PHIST_FORCE_32BIT_GIDX_DEFAULT OFF)
endif()
