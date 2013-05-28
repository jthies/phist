
// type specifier
#ifdef _TP_
#undef _TP_
#endif

#define _TP_ S

// scalar type
#ifdef _ST_
#undef _ST_
#endif

#define _ST_ float

// magnitude (or member) type
#ifdef _MT_
#undef _MT_
#endif
#define _MT_ float

#ifdef _IS_COMPLEX_
#undef _IS_COMPLEX_
#endif

#ifdef _Complex_I
#undef _Complex_I
#endif

#define _Complex_I 0.0f

#ifdef _ZERO_
#undef _ZERO_
#endif

#define _ZERO_ 0.0f

#ifdef _ONE_
#undef _ONE_
#endif

#define _ONE_ 1.0f

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#ifdef _PREF_
#undef _PREF_
#endif

#define _PREF_(name) S ## name

// how to build up the name of a subroutine (void function)
#ifdef _SUBR_
#undef _SUBR_
#endif

#define _SUBR_(name) phist_S ## name

// how to build up the name of a type
#ifdef _TYPE_
#undef _TYPE_
#endif

#define _TYPE_(name) S ## name ## _t

// how to build up the name of a test with two 'template' params
#ifdef _TESTNAME2_
#undef _TESTNAME2_
#endif

#define _TESTNAME2_(name,p1,p2) S ## name ## _ ## p1 ## _ ## p2

// define which gtest macro should be used for comparing floating point numbers
#ifdef ASSERT_REAL_EQ
#undef ASSERT_REAL_EQ
#endif

#define ASSERT_REAL_EQ(expected,actual) ASSERT_FLOAT_EQ(expected,actual)

