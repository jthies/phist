#include "phist_gen_pre.h"

// this may be useful for e.g. generating
// the fortran interfaces or something, to
// avoid littering the preprocessed files
// with lots of system stuff
#ifndef NO_INCLUDES_IN_HEADERS
// let the kernel lib define the complex data type
#include "phist_typedefs.h"
#endif

#define _IS_COMPLEX_
#define _ST_ s_complex_t

// C++ users should use the class phist::ScalarTraits and
// phist_std_typedefs.hpp for all of this:
#ifndef __cplusplus

// type specifier
#define _TP_ 'C'

// scalar type
#define _CMPLX_I_ (0.0f+1.0f*I)
#define _ZERO_ (0.0f+0.0f*I)
#define _ONE_ (1.0f+0.0f*I)
#define _SQRT_(X) csqrtf(X)
#define _CONJ_(X) conjf(X)
#define _ABS_(X) cabsf(X)
#define _REAL_(X) crealf(X)
#define _IMAG_(X) cimagf(X)

#endif


// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) C ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_C ## name

// how to build up the name of a type
#define _TYPE_(name) C ## name ## _t

#include "phist_gen_post.h"
