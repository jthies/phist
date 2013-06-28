#ifndef _IS_COMPLEX_
#define _CONJ_(x) x
#define _REAL_(x) x
#define _IMAG_(x) _ZERO_
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

#ifndef _TP_
#error "macro _TP_ not defined"
#endif

#ifndef _ST_
#error "macro _ST_ not defined"
#endif

#ifndef _CMPLX_I_
#warning "macro _CMPLX_I_ not defined"
#endif

#ifndef _MT_
#error "macro _MT_ not defined"
#endif

#ifndef _CONJ_
#error "macro _CONJ_ not defined"
#endif

#ifndef _SQRT_
#error "macro _SQRT_ not defined"
#endif

#ifndef _ABS_
#error "macro _ABS_ not defined"
#endif

#ifndef _REAL_
#error "macro _REAL_ not defined"
#endif

#ifndef _IMAG_
#error "macro _IMAG_ not defined"
#endif

#ifndef _ZERO_
#error "macro _ZERO_ not defined"
#endif

#ifndef _ONE_
#error "macro _ONE_ not defined"
#endif

#ifndef _PREF_
#error "macro _PREF_ not defined"
#endif

#ifndef _SUBR_
#error "macro _SUBR_ not defined"
#endif

#ifndef _TYPE_
#error "macro _TYPE_ not defined"
#endif

