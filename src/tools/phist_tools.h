/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_TOOLS_H
#define PHIST_TOOLS_H
//! \file phist_tools.h 
#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
#include "phist_ScalarTraits.hpp"
#endif

#ifdef __cplusplus
#include <cstdio>
#include <iostream>
#else
# include <stdio.h>
#endif

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#ifdef PHIST_HAVE_GHOST
#include <ghost.h>
#endif

#endif //DOXYGEN

//! \defgroup tools tools: Utility functions and macros
//!@{

# include "phist_macros.h"

#ifdef __cplusplus

//! \brief get the size of a cache line on the current MPI process,
//! or 1 if the architecture doesn't have a cache.
//!
//! The unit is sizeof(T), so on a typical CPU with 64 bit cache lines
//! phist_cacheline_size<double>()==8 and phist_cacheline_size<char>()==64.
template<typename T>
unsigned int phist_cacheline_size()
{
  static unsigned int clsize=0;
  if (clsize>0) return clsize;
#ifdef PHIST_HAVE_GHOST
  ghost_machine_cacheline_size(&clsize);
# ifdef PHIST_KERNEL_LIB_GHOST
  ghost_type type;
  ghost_type_get(&type);
  if (type==GHOST_TYPE_CUDA) clsize=sizeof(T);
# endif  
  
#else
  PHIST_SOUT(PHIST_WARNING,"guessing cacheline length (64 byte) because ghost is not available.\n");
  clsize=64;
#endif
  clsize/=sizeof(T);
  return clsize;
}

  //! \brief redirect all subsequent output in phist to this output stream
  //! (if this function is not called, the default stream is std:cout)
  //!
  //! Any subsequent call to phist_set_CXX_output_stream or phist_set_C_output_stream
  //! overrides the previous setting.
  //!
  //! All output can be suppressed by calling this function with a nullptr argument.
  void phist_set_CXX_output_stream(std::ostream& ostr);

  //! \brief retrieve the C++ stream to which phistj is printing.
  //!
  //! May be nullptr if a C output stream was set using phist_set_C_output_stream.
  std::ostream* phist_get_CXX_output_stream();

extern "C" {
#endif

#ifdef PHIST_HAVE_MPI
  //! \brief set the standard MPI communicator used to create all subsequent objects
  //! (the one wrapped and returned by phist_comm_create)
  //!
  //! \warning for many kernel libaries (e.g. builtin, ghost, tpetra) we implement 
  //!         a pinning strategy (i.e. bind MPI processes and threads to physical 
  //!         cores of the machine). This is typically done in phist_kernels_init,
  //!         so we recommend setting the communicator *before* that function. If 
  //!         you want to be able to spawn new MPI processes, you may want to set 
  //!         PHIST_TRY_TO_PIN_THREADS=OFF using cmake and use this function to   
  //!         update the communicator.
  void phist_set_default_comm(MPI_Comm new_comm);
  
  //! return the MPI communicator used internally by default
  MPI_Comm phist_get_default_comm();
  
  //! this function should be called from Fortran to set the MPI communicator
  //! analogously to phist_set_default_comm.
  void phist_set_default_comm_f(MPI_Fint f_comm);

  //! \brief this function should be called from Fortran to get a valid MPI communicator
  //! analogously to phist_get_default_comm.
  //!
  //! For a fortran binding of this       
  //! function, see src/kernels/builtin/env_module.f90.
  MPI_Fint phist_get_default_comm_f();
#endif
  //! \brief redirect all subsequent output in phist to this C output stream
  //! (if this function is not called, the default stream is stdout)
  //!
  //! Any subsequent call to phist_set_CXX_output_stream or phist_set_C_output_stream
  //! overrides the previous setting.
  //!
  //! All output can be suppressed by calling this function with a NULL argument.
  void phist_set_C_output_stream(FILE* fp);

  //! return the version of phist as a string (formatted e.g. as "1.7.2")
  const char* phist_version();
  //! return the exact git revision of the phist installation
  const char* phist_git_revision();
  //! \brief return some other useful information on the phist installation as a multi-line 
  //! string, e.g. compilers, compile flags etc.
  const char* phist_install_info();

  //! \brief Get the default output stream (C-style) to which phist is printing.
  //!
  //! May be NULL if a C++ stream was provided using phist_set_CXX_output_stream().
  FILE* phist_get_C_output_stream();

  //! \brief return the name of the library that provides the kernels for this phist installation as a string,
  //! e.g. "builtin","ghost","tpetra".
  const char* phist_kernel_lib();
  
  //! return 1 if mvecs are stored in row-major order, 0 otherwise
  int mvecs_rowmajor();

  //! convert a standard return code from phist into a human-readable string.
  const char* phist_retcode2str(int code);

#ifdef PHIST_HAVE_GHOST
#include <ghost/config.h>
#include <ghost/types.h>
//! for internal use only: convert a ghost error code into a human-readable string
const char* phist_ghost_error2str(ghost_error code);
#endif

# ifdef PHIST_HAVE_MPI

//! \brief pretty-print process-local strings in order.
//!
//! This function should
//! not be used directly but via the wrapper macro PHIST_ORDERED_OUT(...)
int phist_ordered_fprintf(MPI_Comm comm, const char* fmt, ...);

# endif

//!@}

#ifdef __cplusplus
}
#endif

#include "phist_tasks.h"

#endif
