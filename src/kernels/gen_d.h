
// how to build up the name of a subroutine (void function)
#ifdef _SUBROUTINE_
#undef _SUBROUTINE_
#endif

#define _SUBROUTINE_(name) void D ## name

// how to build up the name of a type
#ifdef _TYPENAME_
#undef _TYPENAME_
#endif

#define _TYPENAME_(name) D ## name ## _t

// scalar type
#ifdef _ST_
#undef _ST_
#endif
#define _ST_ double

// magnitude (or member) type
#ifdef _MT_
#undef _MT_
#endif
#define _MT_ double


#ifdef _IS_COMPLEX
#undef _IS_COMPLEX_
#endif


