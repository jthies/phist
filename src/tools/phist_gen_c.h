#include "phist_gen_pre.h"

#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif

#define _IS_COMPLEX_

// type specifier
#define _TP_ C

// scalar type
#ifdef __cplusplus
#define _ST_ std::complex<float>
#define _CMPLX_I_ _ST_(0.0f,1.0f)
#define _ZERO_ _ST_(0.0f,0.0f)
#define _ONE_ _ST_(1.0f,0.0f)
#else
#define _ST_ float complex
#define _CMPLX_I_ (0.0f+1.0f*I)
#define _ZERO_ (0.0f+0.0f*I)
#define _ONE_ (1.0f+0.0f*I)
#endif


// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) C ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_C ## name

// how to build up the name of a type
#define _TYPE_(name) C ## name ## _t

#include "phist_gen_post.h"
