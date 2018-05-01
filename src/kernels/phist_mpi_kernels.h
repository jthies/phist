/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_MPI_KERNELS_H
#define PHIST_MPI_KERNELS_H

#include "phist_config.h"

#ifndef DOXYGEN

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
typedef int MPI_Request;
typedef int MPI_Status;
#endif

/* note: the phist_typedefs.h file is provided in the subdirectory
 where the interface is implemented (e.g. ghost/, tpetra/).
 */
#include "phist_typedefs.h"

#endif /* doxygen */

#ifdef __cplusplus
extern "C" {
#endif

//! \defgroup mpi_kernels communication operations implemented for all kernel libraries in one place
//! \ingroup kernels
//!@{

#define PHIST_CLASSFILE_DEF "phist_mpi_kernels_decl.h"
#include "phist_gen_all.h"

//!@}

#ifdef __cplusplus
} /*extern "C"*/
#endif

#endif
