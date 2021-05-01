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
#include <complex>
#include <cstddef>
#else
//#include <inttypes.h>
#include <complex.h>
#include <stddef.h>
#endif
#endif /* DOXYGEN */

#ifndef DOXYGEN
#include "TpetraCore_config.h"
#else
/* we need something, so we guess the default */
# ifdef PHIST_FORCE_32BIT_GIDX
# define HAVE_TPETRA_INST_INT_LONG
# else
# define HAVE_TPETRA_INST_INT_LONG_LONG
# endif
#endif

#ifdef __cplusplus
//! single precision complex type
using phist_s_complex = std::complex<float>;
//! double precision complex type
using phist_d_complex = std::complex<double>;
#else
typedef float complex phist_s_complex;
typedef double complex phist_d_complex;
#endif

//! type of global indices
#ifdef PHIST_FORCE_32BIT_GIDX
#if defined(HAVE_TPETRA_INST_INT_INT)
typedef int phist_gidx;
#elif defined(HAVE_TPETRA_INST_INT_LONG)
typedef long phist_gidx;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED)
typedef unsigned phist_gidx;
#elif defined(HAVE_TPETRA_INST_INT_UNSIGNED_LONG)
typedef unsigned long phist_gidx;
#else
#error "You requested 32-bit global indices, but the Tpetra installation does not instantiate any int/long global index types supported here."
#endif
#define PRgidx "d"
#else
#ifdef HAVE_TPETRA_INST_INT_LONG_LONG
typedef long long phist_gidx;
#else
#error "You requested 64-bit global indices, but the Tpetra installation does not instantiate any long_long global index types supported here."
#endif
#define PRgidx "lld"
#endif

// we want ptrdiff_t (aka long long int on 64 bit systems) as local index,
// but a bug in Trilinos prevents us from using it right now. until then,
// we use int as local index type

//! type of node-local indices
typedef int phist_lidx;
#define PRlidx "d"

#ifndef DOXYGEN
#include "phist_void_aliases.h"
#endif //DOXYGEN

#endif
