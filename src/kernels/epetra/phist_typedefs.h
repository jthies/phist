#ifndef PHIST_TYPEDEFS_H
#define PHIST_TYPEDEFS_H

#ifdef __cplusplus

#ifndef NO_INCLUDES_IN_HEADERS
//#include <cinttypes>
#include <complex>
//! single precision complex type
typedef std::complex<float> s_complex_t;
typedef std::complex<double> d_complex_t;
#endif

#else

#ifndef NO_INCLUDES_IN_HEADERS
//#include <inttypes.h>
#include <complex.h>
typedef float complex s_complex_t;
typedef double complex d_complex_t;
#endif

#endif


//! type of node-local indices
typedef int lidx_t;

//! type of global indices
typedef int gidx_t;

//#define PRlidx=PRId32
//#define PRgidx=PRId32
#define PRlidx "d"
#define PRgidx "d"

#include "phist_void_aliases.h"

#endif
