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
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <inttypes.h>

//! complex data types
#ifdef __cplusplus
#include <complex>
typedef std::complex<float> phist_s_complex;
typedef std::complex<double> phist_d_complex;
#else
#include <complex.h>
typedef  complex float phist_s_complex;
typedef  complex double phist_d_complex;
#endif
#endif /* DOXYGEN */
//! type of node-local indices
typedef int32_t phist_lidx;

//! type of global indices
typedef int64_t phist_gidx;

#define PRlidx "d"
#define PRgidx "ld"

#ifdef PHIST_FORCE_INT_GIDX
# warning "neglecting config option PHIST_FORCE_INT_GIDX with builtin kernels!"
#endif

#include "phist_void_aliases.h"

#endif
