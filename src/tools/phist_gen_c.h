#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif

// type specifier
#ifdef _TP_
#undef _TP_
#endif

#define _TP_ C

// scalar type
#ifdef _ST_
#undef _ST_
#endif

#ifdef __cplusplus
#define _ST_ std::complex<float>
#ifndef _Complex_I
#define _Complex_I _ST_(0.0,1.0)
#endif
#else
#define _ST_ _Complex float
#endif


// magnitude (or member) type
#ifdef _MT_
#undef _MT_
#endif
#define _MT_ float

#ifdef _IS_COMPLEX_
#undef _IS_COMPLEX_
#endif

#define _IS_COMPLEX_


// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#ifdef _PREF_
#undef _PREF_
#endif

#define _PREF_(name) C ## name

// how to build up the name of a subroutine (void function)
#ifdef _SUBR_
#undef _SUBR_
#endif

#define _SUBR_(name) phist_C ## name

// how to build up the name of a type
#ifdef _TYPE_
#undef _TYPE_
#endif

#define _TYPE_(name) C ## name ## _t

