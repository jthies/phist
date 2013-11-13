#include "phist_gen_pre.h"


#define IS_COMPLEX
#define _ST_ s_complex_t

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus

// type specifier
#define _TP_ 'C'

#define _CMPLX_I_ (0.0f+1.0f*I)
#define ZERO (0.0f+0.0f*I)
#define ONE (1.0f+0.0f*I)
#define SQRT(X) csqrtf(X)
#define MSQRT(X) sqrtf(X)
#define CONJ(X) conjf(X)
#define ABS(X) cabsf(X)
#define MABS(X) absf(X)
#define REAL(X) crealf(X)
#define IMAG(X) cimagf(X)

#endif


// adds type prefix to a specifier, e.g. PREFIX(gemm) -> Dgemm
#define PREFIX(name) C ## name

// how to build up the name of a subroutine (void function)
#define SUBR(name) phist_C ## name

// how to build up the name of a type
#define TYPE(name) C ## name ## _t

#include "phist_gen_post.h"
