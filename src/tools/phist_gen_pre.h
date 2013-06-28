// this file undefs all the macros defined by the phist_gen_* 
// headers

// type specifier
#ifdef _TP_
#undef _TP_
#endif

// scalar type
#ifdef _ST_
#undef _ST_
#endif

// imaginary unit
#ifdef _CMPLX_I_
#undef _CMPLX_I_
#endif

// magnitude (or member) type
#ifdef _MT_
#undef _MT_
#endif

// how to take the conjugate of a scalar
#ifdef _CONJ_
#undef _CONJ_
#endif

// squareroot
#ifdef _SQRT_
#undef _SQRT_
#endif

// wether _ST_ is a real or complex type
#ifdef _IS_COMPLEX_
#undef _IS_COMPLEX_
#endif

// wether _ST_ is a single or double precision type
#ifdef _IS_DOUBLE_
#undef _IS_DOUBLE_
#endif

// scalar 0 in the type _ST_
#ifdef _ZERO_
#undef _ZERO_
#endif

// scalar 1 in the type _ST_
#ifdef _ONE_
#undef _ONE_
#endif

// adds type prefix to a specifier, e.g. _PREF_(gemm) -> Dgemm
#ifdef _PREF_
#undef _PREF_
#endif

// how to build up the name of a subroutine (void function)
#ifdef _SUBR_
#undef _SUBR_
#endif

// how to build up the name of a type
#ifdef _TYPE_
#undef _TYPE_
#endif

// macro used in googletest assertions
#ifdef ASSERT_REAL_EQ
#undef ASSERT_REAL_EQ
#endif
