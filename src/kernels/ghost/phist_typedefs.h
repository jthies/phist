#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "phist_config.h"
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
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
#else
//#include <inttypes.h>
#include <complex.h>
typedef  complex float s_complex_t;
typedef  complex double d_complex_t;
#endif

#include "ghost/config.h"
#include "ghost/types.h"

//! type of node-local indices. 
typedef ghost_lidx_t lidx_t;

//! type of global indices
typedef ghost_gidx_t gidx_t;

#ifdef GHOST_HAVE_LONGIDX_LOCAL
#define PRlidx "lld"
#else
#define PRlidx "d"
#endif

#ifdef GHOST_HAVE_LONGIDX_GLOBAL
#define PRgidx "lld"
#else
#define PRgidx "d"
#endif

#include "phist_void_aliases.h"

#endif
