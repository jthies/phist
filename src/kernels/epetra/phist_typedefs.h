#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#ifdef __cplusplus
#ifndef NO_INCLUDES_IN_HEADERS
#include <complex>
//! single precision complex type
typedef std::complex<float> s_complex_t;
#else
#include <complex.h>
typedef float complex s_complex_t;
typedef double complex d_complex_t;
#endif
#endif


//! type of node-local indices
typedef int lidx_t;

//! type of global indices
typedef int gidx_t;

#include "phist_void_aliases.h"

#endif
