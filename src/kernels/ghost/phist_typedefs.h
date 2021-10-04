/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
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

#include "phist_macros.h"
//! complex data types
#ifdef __cplusplus
//#include <cinttypes>
#include <complex>
//!  TODO - ghost has its own
//! C++ template, but I'm not sure they are usable
//! with Belos and TSQR. The data layout should
//! be the same, however.
typedef std::complex<float> phist_s_complex;
typedef std::complex<double> phist_d_complex;
#else
//#include <inttypes.h>
#include <complex.h>
typedef  complex float phist_s_complex;
typedef  complex double phist_d_complex;
#endif

#include <ghost/config.h>
#include <ghost/types.h>

#endif /* DOXYGEN */

//! type of node-local indices. 
typedef ghost_lidx phist_lidx;

//! type of global indices
typedef ghost_gidx phist_gidx;

#ifdef GHOST_IDX64_LOCAL
#define PRlidx "lld"
#else
#define PRlidx "d"
#endif

#ifdef GHOST_IDX64_GLOBAL
#define PRgidx "lld"
#else
#define PRgidx "d"
#endif

#ifdef PHIST_FORCE_32BIT_GIDX
# warning "neglecting config option PHIST_FORCE_32BIT_GIDX with ghost kernels, choice depends on GHOOST installation!"
#endif

#ifndef DOXYGEN
#include "phist_void_aliases.h"
#endif //DOXYGEN

#endif
