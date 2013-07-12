#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#ifndef NO_INCLUDES_IN_HEADERS
#ifdef __cplusplus
#include <complex>
#include <cstddef>
#else
#include <complex.h>
#include <stddef.h>
#endif
#endif

#ifdef __cplusplus
//! single precision complex type
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
//! type of global indices
typedef std::ptrdiff_t gidx_t;
#else
typedef float complex s_complex_t;
typedef double complex d_complex_t;
//! type of global indices
typedef ptrdiff_t gidx_t;
#endif


// we want ptrdiff_t (aka long long int on 64 bit systems) as local index,
// but a bug in Trilinos prevents us from using it right now. until then,
// we use int as local index type

//! type of node-local indices
typedef int lidx_t;

#include "phist_void_aliases.h"

#endif
