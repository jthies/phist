#include "phist_gen_pre.h"

// this may be useful for e.g. generating
// the fortran interfaces or something, to
// avoid littering the preprocessed files
// with lots of system stuff
#ifndef NO_INCLUDES_IN_HEADERS
#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif
#endif

// type specifier
#define _TP_ 'Z'

#define _IS_DOUBLE_
#define _IS_COMPLEX_

// scalar type
#ifdef __cplusplus
#define _ST_ std::complex<double>
#define _CMPLX_I_ (_ST_(0.0,1.0))
#define _ZERO_ (_ST_(0.0,0.0))
#define _ONE_ (_ST_(1.0,0.0))
#define _SQRT_(X) std::sqrt(X)
#define _ABS_(X) std::abs(X)
#define _CONJ_(X) std::conj(X)
#define _REAL_(X) std::real(X)
#define _IMAG_(X) std::imag(X)
#else
#define _ST_ double complex
#define _CMPLX_I_ (0.0+1.0*I)
#define _ZERO_ (0.0+0.0*I)
#define _ONE_ (1.0+0.0*I)
#define _SQRT_(X) csqrt(X)
#define _ABS_(X) cabs(X)
#define _CONJ_(X) conj(X)
#define _REAL_(X) creal(X)
#define _IMAG_(X) cimag(X)
#endif

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) Z ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_Z ## name

// how to build up the name of a type
#define _TYPE_(name) Z ## name ## _t

#include "phist_gen_post.h"
