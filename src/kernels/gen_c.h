#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif
// how to build up the name of a subroutine (void function)
#ifdef _SUBROUTINE_
#undef _SUBROUTINE_
#endif

#define _SUBROUTINE_(name) void essex_C ## name

// how to build up the name of a type
#ifdef _TYPE_
#undef _TYPE_
#endif

#define _TYPE_(name) C ## name ## _t

// scalar type
#ifdef _ST_
#undef _ST_
#endif

#ifdef __cplusplus
#define _ST_ std::complex<float>
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

