#include "phist_gen_pre.h"

#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif

// type specifier
#define _TP_ Z

#define _IS_DOUBLE_
#define _IS_COMPLEX_

// scalar type
#ifdef __cplusplus
#define _ST_ std::complex<double>
#define _Complex_I (_ST_(0.0,1.0))
#else
#define _ST_ _Complex double
#endif

#define _ZERO_ _ST_(0.0,0.0)

#define _ONE_ _ST_(1.0,0.0)

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#define _PREF_(name) Z ## name

// how to build up the name of a subroutine (void function)
#define _SUBR_(name) phist_Z ## name

// how to build up the name of a type
#define _TYPE_(name) Z ## name ## _t

#include "phist_gen_post.h"
