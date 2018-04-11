/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"

#ifndef DOXYGEN
/* needs to be included before system headers for some intel compilers+mpi */
# ifdef PHIST_HAVE_MPI
# include <mpi.h>
# endif
#endif
#ifdef __cplusplus

#ifndef DOXYGEN
//#include <cinttypes>
#include <complex>
//! single precision complex type
typedef std::complex<float> phist_s_complex;
typedef std::complex<double> phist_d_complex;
#endif

#else

#ifndef DOXYGEN
//#include <inttypes.h>
#include <complex.h>
typedef float complex phist_s_complex;
typedef double complex phist_d_complex;
#endif

#endif

#include "Epetra_config.h"

//! type of node-local indices
typedef int phist_lidx;

//! type of global indices
#if defined(EPETRA_NO_64BIT_GLOBAL_INDICES)||defined(PHIST_FORCE_32BIT_GIDX)
typedef int phist_gidx;
#define PRgidx "ld"
#else
typedef long long int phist_gidx;
#define PRgidx "d"
#endif

#define PRlidx "d"

#include "phist_void_aliases.h"

#endif
