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
#ifdef _Complex_I
#undef _Complex_I
#endif
#define _Complex_I _ST_(0.0f,1.0f)
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

#ifdef _ZERO_
#undef _ZERO_
#endif

#define _ZERO_ _ST_(0.0f,0.0f)

#ifdef _ONE_
#undef _ONE_
#endif

#define _ONE_ _ST_(1.0f,0.0f)

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


// how to build up the name of a test with two 'template' params
#ifdef _TESTNAME2_
#undef _TESTNAME2_
#endif

#define _TESTNAME2_(name,p1,p2) C ## name ## _ ## p1 ## _ ## p2

#include "phist_macros.h"
