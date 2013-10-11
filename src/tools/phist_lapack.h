#ifndef PHIST_LAPACK_H
#define PHIST_LAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Slap_cmplx_t {
float re;
float im;
} Slap_cmplx_t;

typedef struct Dlap_cmplx_t {
double re;
double im;
} Dlap_cmplx_t;

// TODO - cmake/blas/lapack integration.
// this is a platform dependent macro, we should have CMake determine
// how to define the name of a fortran 77 routine
#define LAPACK_SUBR(NAME,name) name ## _

// we add lapack subroutines as we go along implementing things.

///////////////////////////////////////////////////////////////////////////////////////////
// XSTEQR - QR decomposition of symmetric tridiagonal matrices                           //
///////////////////////////////////////////////////////////////////////////////////////////


// QR decomposition of a real symmetric tridiagonal matrix
#define SSTEQR LAPACK_SUBR(SSTEQR,ssteqr)
void SSTEQR(const char*, const int* n, float* D, float* E, float* Z, const int* ldz, float* work, int* info);
// QR decomposition of a real symmetric tridiagonal matrix
#define DSTEQR LAPACK_SUBR(DSTEQR,dsteqr)
void DSTEQR(const char*, const int* n, double* D, double* E, double* Z, const int* ldz, double* work, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEES - Schur decomposition                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

// Schur decomposition of a real non-symmetric matrix
#define SGEES LAPACK_SUBR(SGEES,sgees)
void DGEES(const char*, const char*, int (*select)(float*, float*), const int* n, 
float* a, const int* lda, int* sdim, float* wr, float* wi,
 float* vs, const int* ldvs, float* work, const int* lwork, int* bwork, int* info);
// Schur decomposition of a real non-symmetric matrix
#define DGEES LAPACK_SUBR(DGEES,dgees)
void DGEES(const char*, const char*, int (*select)(double*, double*), const int* n, 
double* a, const int* lda, int*sdim, double* wr, double* wi,
 double* vs, const int* ldvs, double* work, const int* lwork, int* bwork, int* info);


// Schur decomposition of a complex non-hermitian matrix 
#define CGEES LAPACK_SUBR(CGEES,cgees)
void ZGEES(const char*, const char*, int (*select)(Slap_cmplx_t*), 
const int* n, Slap_cmplx_t* a, const int* lda, int* sdim,
Slap_cmplx_t* w, Slap_cmplx_t* vs, const int* ldvs, Slap_cmplx_t* work, 
const int* lwork, float* rwork, int* bwork, int* info);
// Schur decomposition of a complex non-hermitian matrix 
#define ZGEES LAPACK_SUBR(ZGEES,zgees)
void ZGEES(const char*, const char*, int (*select)(Dlap_cmplx_t*), 
const int* n, Dlap_cmplx_t* a, const int* lda, int* sdim,
Dlap_cmplx_t* w, Dlap_cmplx_t* vs, const int* ldvs, Dlap_cmplx_t* work, 
const int* lwork, double* rwork, int* bwork, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRSEN - reorder Schur form                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

//! reorder real Schur decomposition
#define STRSEN LAPACK_SUBR(STRSEN,strsen)
void STRSEN(const char* job, const char* jobvs, const int *select, const int *n, float *t, 
const int *ldt, float *q, const int *ldq, float *wr, float *wi, int *m, float *s, float *sep, 
float *work, const int *lwork, int *iwork, const int *liwork, int *info);

//! reorder real Schur decomposition
#define DTRSEN LAPACK_SUBR(DTRSEN,dtrsen)
void DTRSEN(const char* job, const char* jobvs, const int *select, const int *n, double *t, 
const int *ldt, double *q, const int *ldq, double *wr, double *wi, int *m, double *s, double *sep, 
double *work, const int *lwork, int *iwork, const int *liwork, int *info);

//! reorder complex Schur decomposition
#define CTRSEN LAPACK_SUBR(CTRSEN,ctrsen)
void CTRSEN(const char* job, const char* jobvs, const int *select, const int *n, 
Slap_cmplx_t *T, const int *ldt, Slap_cmplx_t *Q, const int *ldq, Slap_cmplx_t *W, 
const int *m, double *S, double *sep, Slap_cmplx_t *work, const int *lwork, int *info);

//! reorder complex Schur decomposition
#define ZTRSEN LAPACK_SUBR(ZTRSEN,ztrsen)
void ZTRSEN(const char* job, const char* jobvs, const int *select, const int *n, 
Dlap_cmplx_t *T, const int *ldt, Dlap_cmplx_t *Q, const int *ldq, Dlap_cmplx_t *W, 
const int *m, double *S, double *sep, Dlap_cmplx_t *work, const int *lwork, int *info);

#ifdef __cplusplus
}
#endif

#endif
