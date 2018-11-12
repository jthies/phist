include(LibFindMacros)

libfind_pkg_check_modules(PETSC PETSc)
libfind_process(PETSC)

# Figure out which data type is supported by the PETSc installation found.
# Note that a given PETSc installation supports exactly one data type (S/C/D/Z)
# at a time.
include(CheckSymbolExists)
check_symbol_exists(PETSC_USE_REAL_DOUBLE "petscsys.h" PETSC_USE_REAL_DOUBLE)
check_symbol_exists(PETSC_USE_REAL_SINGLE "petscsys.h" PETSC_USE_REAL_SINGLE)
check_symbol_exists(PETSC_USE_COMPLEX "petscsys.h" PETSC_USE_COMPLEX)

if (PETSC_USE_REAL_SINGLE AND PETSC_USE_REAL_DOUBLE)
  message(WARNING "Found PETSc installation, but it unexpectedly defines both "
                  "'PETSC_USE_REAL_SINGLE' and 'PETSC_USE_REAL_DOUBLE'."
                  "Will continue with both data types enabled.")
endif()
