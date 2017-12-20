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

#ifdef __cplusplus
//#include <cinttypes>
//#include <complex>
#include <Kokkos_Complex.hpp>
#include <cstddef>
#else
//#include <inttypes.h>
#include <complex.h>
#include <stddef.h>
#endif
#endif /* DOXYGEN */

#ifdef __cplusplus
//! single precision complex type
using phist_s_complex = std::complex<float>;
//! double precision complex type
using phist_d_complex = std::complex<double>;
//! type of global indices
using phist_gidx = std::ptrdiff_t;
#else
typedef float complex phist_s_complex;
typedef double complex phist_d_complex;
//! type of global indices
#ifdef PHIST_FORCE_INT_GIDX
typedef int phist_gidx;
#define PRgidx "d"
#else
typedef ptrdiff_t phist_gidx;
#define PRgidx "ld"
#endif
#endif

// we want ptrdiff_t (aka long long int on 64 bit systems) as local index,
// but a bug in Trilinos prevents us from using it right now. until then,
// we use int as local index type

//! type of node-local indices
typedef int phist_lidx;
#define PRlidx "d"

#include "phist_void_aliases.h"

#endif
