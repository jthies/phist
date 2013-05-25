
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

#ifdef _IS_COMPLEX
#undef _IS_COMPLEX_
#endif


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

