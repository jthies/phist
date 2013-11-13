#include "phist_gen_pre.h"


// type specifier
#define _TP_ 'Z'
#define _ST_ d_complex_t

#define _IS_DOUBLE_
#define IS_COMPLEX

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus
#define _CMPLX_I_ (0.0+1.0*I)
#define ZERO (0.0+0.0*I)
#define ONE (1.0+0.0*I)
#define SQRT(X) csqrt(X)
#define MSQRT(X) sqrt(X)
#define ABS(X) cabs(X)
#define MABS(X) abs(X)
#define CONJ(X) conj(X)
#define REAL(X) creal(X)
#define IMAG(X) cimag(X)
#endif

// adds type prefix to a specifier, e.g. PREFIX(gemm) -> Dgemm
#define PREFIX(name) Z ## name

// how to build up the name of a subroutine (void function)
#define SUBR(name) phist_Z ## name

// how to build up the name of a type
#define TYPE(name) Z ## name ## _t

#include "phist_gen_post.h"
