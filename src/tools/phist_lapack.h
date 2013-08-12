#ifndef PHIST_LAPACK_H
#define PHIST_LAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

// this is a platform dependent macro, we should have CMake determine
// how to define the name of a fortran 77 routine
#define _LAPACKSUBR(NAME,name) name ## _

// we add lapack subroutines as we go along implementing things.

// QR decomposition of a real symmetric tridiagonal matrix
#define SSTEQR _LAPACKSUBR(SSTEQR,ssteqr)
void SSTEQR(const char*, const int* n, float* D, float* E, float* Z, const int* ldz, float* work, int* info);
// QR decomposition of a real symmetric tridiagonal matrix
#define DSTEQR _LAPACKSUBR(DSTEQR,dsteqr)
void DSTEQR(const char*, const int* n, double* D, double* E, double* Z, const int* ldz, double* work, int* info);

#ifdef __cplusplus
}
#endif

#endif
