# Figure out which data types are supported by the Trilinos installation found.
include(CheckSymbolExists)

set(CMAKE_REQUIRED_INCLUDES ${Teuchos_INCLUDE_DIRS} ${Tpetra_INCLUDE_DIRS})
check_symbol_exists(HAVE_TEUCHOS_COMPLEX "Teuchos_config.h" HAVE_TEUCHOS_COMPLEX)
check_symbol_exists(HAVE_TPETRA_INST_DOUBLE "TpetraCore_config.h" HAVE_TPETRA_INST_DOUBLE)
check_symbol_exists(HAVE_TPETRA_INST_FLOAT "TpetraCore_config.h" HAVE_TPETRA_INST_FLOAT)
check_symbol_exists(HAVE_TPETRA_INST_COMPLEX_DOUBLE "TpetraCore_config.h" HAVE_TPETRA_INST_COMPLEX_DOUBLE)
check_symbol_exists(HAVE_TPETRA_INST_COMPLEX_FLOAT "TpetraCore_config.h" HAVE_TPETRA_INST_COMPLEX_FLOAT)

if (PHIST_KERNEL_LIB STREQUAL "tpetra")
  if (NOT HAVE_TPETRA_INST_DOUBLE)
    message(FATAL_ERROR "PHIST can't be compiled with Tpetra if Tpetra does not support real double data "
                        "(TPETRA_INST_DOUBLE is not defined)")
  endif()
  
  set(PHIST_SUPPORT_SP OFF)
  set(PHIST_SUPPORT_COMPLEX OFF)
  
  if (HAVE_TPETRA_INST_COMPLEX_DOUBLE)
    set(PHIST_SUPPORT_COMPLEX ON)
  endif()
  
  if (HAVE_TPETRA_INST_FLOAT)
    set(PHIST_SUPPORT_SP ON)
  endif()
  
  if (PHIST_SUPPORT_COMPLEX AND NOT HAVE_TEUCHOS_COMPLEX)
    message(WARNING "disabling complex arithmetic support (PHIST_SUPPORT_COMPLEX) because Teuchos is compiled"
                    "without HAVE_TEUCHOS_COMPLEX")
    set(PHIST_SUPPORT_COMPLEX OFF)
  endif()

  if (PHIST_SUPPORT_COMPLEX AND PHIST_SUPPORT_SP AND NOT HAVE_TPETRA_INST_COMPLEX_FLOAT)
    message(WARNING "disabling single precision support (PHIST_SUPPORT_SP) because Tpetra is compiled with"
                    "HAVE_TPETRA_INST_COMPLEX_DOUBLE and HAVE_TPETRA_INST_FLOAT, but without HAVE_TPETRA_INST_COMPLEX_FLOAT")
    set(PHIST_SUPPORT_SP OFF)
  endif()
endif()
