# this file tries to figure out a compiler flag to enable host-specific
# code generation and returns it in ARCH_FLAG (e.g. -march=native, -xHOST, ...)
# If the variable PHIST_HOST_OPTIMIZE is not set, ARCH_FLAG is left blank and
# no host-specific code (like SIMD operations) will be generated.

include(CheckCCompilerFlag)
include(CheckFortranCompilerFlag)
include(CheckCXXCompilerFlag)

cmake_push_check_state() # Save variables

function(CHECK_COMPILER_FLAG lang flag var)
  if (lang STREQUAL "C")
    CHECK_C_COMPILER_FLAG(${flag} ${var})
  elseif (lang STREQUAL "CXX")
    CHECK_CXX_COMPILER_FLAG(${flag} ${var})
  elseif (lang STREQUAL "Fortran")
    CHECK_Fortran_COMPILER_FLAG(${flag} ${var})
  endif()
endfunction()

if(PHIST_HOST_OPTIMIZE)

  # typical flags for this (at least for GNU, Intel, CLANG):
  set(FLAG_LIST -march=native -mcpu=native -xHOST)

  foreach (LANG C CXX Fortran)
    foreach (FLAG ${FLAG_LIST})
      CHECK_COMPILER_FLAG(${LANG} ${FLAG} ${LANG}_ARCH_FLAG_FOUND)
      if (${LANG}_ARCH_FLAG_FOUND)
        set(${LANG}_ARCH_FLAG ${FLAG} CACHE STRING "${LANG} compiler flag to enable host-specific code generation")
        break()
      endif()
      list(APPEND CMAKE_${LANG}_FLAGS ${${LANG}_ARCH_FLAG})
    endforeach()
  endforeach()
endif()

cmake_pop_check_state() # Recover variables
