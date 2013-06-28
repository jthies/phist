
#ifdef _LAPACK_SUBR_
#undef _LAPACK_SUBR_
#endif

#ifndef _IS_DOUBLE_
#ifndef _IS_COMPLEX_
#define _LAPACK_SUBR_(NAME,name) s ## name ## _
#else
#define _LAPACK_SUBR_(NAME,name) c ## name ## _
#endif
#else
#ifndef _IS_COMPLEX_
#define _LAPACK_SUBR_(NAME,name) d ## name ## _
#else
#define _LAPACK_SUBR_(NAME,name) z ## name ## _
#endif
#endif

#ifdef STEQR
#undef STEQR
#endif
#define STEQR _LAPACK_SUBR_(STEQR,steqr)

void STEQR(const char*, const int* n, _ST_* D, _ST_* E, _ST_* Z, const int* ldz, _ST_* work, int* info);
