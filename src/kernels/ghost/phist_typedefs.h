#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

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
typedef ghost_idx_t lidx_t;

//! type of global indices
typedef ghost_idx_t gidx_t;

#ifdef GHOST_HAVE_LONGIDX
#define PRlidx "lld"
#define PRgidx "lld"
#else
#define PRlidx "d"
#define PRgidx "d"
#endif

#include "phist_void_aliases.h"

#endif
