#ifdef _IS_COMPLEX_
#define _CONJ_(x) std::conj(x)
#else
#define _CONJ_(x) x
#define _CMPLX_I_ _ZERO_
#endif

#ifdef _IS_DOUBLE_
#define _MT_ double
#define ASSERT_REAL_EQ(expected,actual) ASSERT_DOUBLE_EQ(expected,actual)
#else
#define _MT_ float
#define ASSERT_REAL_EQ(expected,actual) ASSERT_FLOAT_EQ(expected,actual)
#endif

// check wether all necessary macros are there:

// type specifier
#ifndef _TP_
#error "macro _TP_ not defined"
#endif

// scalar type
#ifndef _ST_
#error "macro _ST_ not defined"
#endif

// imaginary unit
#ifndef _CMPLX_I_
#warning "macro _CMPLX_I_ not defined"
#endif

// magnitude (or member) type
#ifndef _MT_
#error "macro _MT_ not defined"
#endif

// how to take the conjugate of a scalar
#ifndef _CONJ_
#error "macro _CONJ_ not defined"
#endif

// scalar 0 in the type _ST_
#ifndef _ZERO_
#error "macro _ZERO_ not defined"
#endif

// scalar 1 in the type _ST_
#ifndef _ONE_
#error "macro _ONE_ not defined"
#endif

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#ifndef _PREF_
#error "macro _PREF_ not defined"
#endif

// how to build up the name of a subroutine (void function)
#ifndef _SUBR_
#error "macro _SUBR_ not defined"
#endif

// how to build up the name of a type
#ifndef _TYPE_
#error "macro _TYPE_ not defined"
#endif

