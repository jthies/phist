#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#include "ghost_types.h"

//! complex data types
#ifdef __cplusplus
#include <complex>
//!  TODO - ghost has its own
//! C++ template, but I'm not sure they are usable
//! with Belos and TSQR. The data layout should
//! be the same, however.
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
#else
#include <complex.h>
typedef  complex float s_complex_t;
typedef  complex double d_complex_t;
#endif

//! type of node-local indices
typedef ghost_vidx_t lidx_t;

//! type of global indices
typedef long long int gidx_t;

#include "phist_void_aliases.h"

#endif
