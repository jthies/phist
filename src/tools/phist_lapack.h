#ifndef PHIST_LAPACK_H
#define PHIST_LAPACK_H

typedef struct Sblas_cmplx_t {
float re;
float im;
} Sblas_cmplx_t;

typedef struct Dblas_cmplx_t {
double re;
double im;
} Dblas_cmplx_t;


#ifdef __cplusplus
extern "C" {
#endif

// TODO - cmake/blas/lapack integration.
// this is a platform dependent macro, we should have CMake determine
// how to define the name of a fortran 77 routine
#define LAPACK_SUBR(NAME,name) name ## _
#define BLAS_SUBR(NAME,name) name ## _

// we add lapack subroutines as we go along implementing things.

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEMM - general dense matrix-matrix multiplicaiton                               //
///////////////////////////////////////////////////////////////////////////////////////////
#define SGEMM BLAS_SUBR(SGEMM,sgemm)
void SGEMM(const char*, const char*, const int*, const int*, const int*,
    const float*, const float*, const int*, const float*, const int*, const float*,
    float *, const int*, int*);
#define DGEMM BLAS_SUBR(DGEMM,dgemm)
void DGEMM(const char*, const char*, const int*, const int*, const int*,
    const double*, const double*, const int*, const double*, const int*, const double*,
    double *, const int*, int*);
#define CGEMM BLAS_SUBR(CGEMM,cgemm)
void CGEMM(const char*, const char*, const int*, const int*, const int*,
    const Sblas_cmplx_t*, const Sblas_cmplx_t*, const int*, const Sblas_cmplx_t*, const int*, const Sblas_cmplx_t*,
    Sblas_cmplx_t *, const int*, int*);
#define ZGEMM BLAS_SUBR(ZGEMM,zgemm)
void ZGEMM(const char*, const char*, const int*, const int*, const int*,
    const Dblas_cmplx_t*, const Dblas_cmplx_t*, const int*, const Dblas_cmplx_t*, const int*, const Dblas_cmplx_t*,
    Dblas_cmplx_t *, const int*, int*);

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
void SGEES(const char*, const char*, int (*select)(float*, float*), const int* n, 
float* a, const int* lda, int* sdim, float* wr, float* wi,
 float* vs, const int* ldvs, float* work, const int* lwork, int* bwork, int* info);
// Schur decomposition of a real non-symmetric matrix
#define DGEES LAPACK_SUBR(DGEES,dgees)
void DGEES(const char*, const char*, int (*select)(double*, double*), const int* n, 
double* a, const int* lda, int*sdim, double* wr, double* wi,
 double* vs, const int* ldvs, double* work, const int* lwork, int* bwork, int* info);


// Schur decomposition of a complex non-hermitian matrix 
#define CGEES LAPACK_SUBR(CGEES,cgees)
void CGEES(const char*, const char*, int (*select)(Sblas_cmplx_t*), 
const int* n, Sblas_cmplx_t* a, const int* lda, int* sdim,
Sblas_cmplx_t* w, Sblas_cmplx_t* vs, const int* ldvs, Sblas_cmplx_t* work, 
const int* lwork, float* rwork, int* bwork, int* info);
// Schur decomposition of a complex non-hermitian matrix 
#define ZGEES LAPACK_SUBR(ZGEES,zgees)
void ZGEES(const char*, const char*, int (*select)(Dblas_cmplx_t*), 
const int* n, Dblas_cmplx_t* a, const int* lda, int* sdim,
Dblas_cmplx_t* w, Dblas_cmplx_t* vs, const int* ldvs, Dblas_cmplx_t* work, 
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
Sblas_cmplx_t *T, const int *ldt, Sblas_cmplx_t *Q, const int *ldq, Sblas_cmplx_t *W, 
const int *m, float *S, float *sep, Sblas_cmplx_t *work, const int *lwork, int *info);

//! reorder complex Schur decomposition
#define ZTRSEN LAPACK_SUBR(ZTRSEN,ztrsen)
void ZTRSEN(const char* job, const char* jobvs, const int *select, const int *n, 
Dblas_cmplx_t *T, const int *ldt, Dblas_cmplx_t *Q, const int *ldq, Dblas_cmplx_t *W, 
const int *m, double *S, double *sep, Dblas_cmplx_t *work, const int *lwork, int *info);

///////////////////////////////////////////////////////////////////////////////////////////
//      TREVC - eigenvalues and -vectors of the Schur form                               //
///////////////////////////////////////////////////////////////////////////////////////////

#define STREVC LAPACK_SUBR(STREVC,strevc)
void STREVC(const char* side, const char* howmny, int* select, const int* n, 
const float* T, const int* ldt, float* vl, const int* ldvl, float* vr, const int* ldvr, 
const int* mm, int* m, float* work, int* info);

#define DTREVC LAPACK_SUBR(DTREVC,dtrevc)
void DTREVC(const char* side, const char* howmny, int* select, const int* n, 
const double* T, const int* ldt, double* vl, const int* ldvl, double* vr, const int* ldvr, 
const int* mm, int* m, double* work, int* info);

#define CTREVC LAPACK_SUBR(CTREVC,ctrevc)
  void CTREVC(const char* side, const char* howmny, int* select, const int* n, 
const Sblas_cmplx_t* t, const int* ldt, Sblas_cmplx_t* vl, const int* ldvl, Sblas_cmplx_t* vr, const int* ldvr, 
const int* mm, int* m, Sblas_cmplx_t* work, float* rwork, int* info);

#define ZTREVC LAPACK_SUBR(ZTREVC,ztrevc)
void ZTREVC(const char* side, const char* howmny, int* select, const int* n, 
const Dblas_cmplx_t* T, const int* ldt, Dblas_cmplx_t* vl, const int* ldvl, Dblas_cmplx_t* vr, const int* ldvr, 
const int* mm, int* m, Dblas_cmplx_t* work, double* rwork, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRS - solve triangular linear system                                          //
///////////////////////////////////////////////////////////////////////////////////////////
#define STRTRS LAPACK_SUBR(STRTRS,strtrs)
void STRTRS(const char* uplo, const char* trans, const char* diag, const int* n, const int* nrhs, 
const float* a, const int* lda, float* b, const int* ldb, int* info);

#define DTRTRS LAPACK_SUBR(DTRTRS,dtrtrs)
void DTRTRS(const char* uplo, const char* trans, const char* diag, const int* n, const int* nrhs, 
const double* a, const int* lda, double* b, const int* ldb, int* info);

#define CTRTRS LAPACK_SUBR(CTRTRS,ctrtrs)
void CTRTRS(const char* uplo, const char* trans, const char* diag, const int* n, const int* nrhs, 
const Sblas_cmplx_t* a, const int* lda, Sblas_cmplx_t* b, const int* ldb, int* info);

#define ZTRTRS LAPACK_SUBR(ZTRTRS,ztrtrs)
void ZTRTRS(const char* uplo, const char* trans, const char* diag, const int* n, const int* nrhs, 
const Dblas_cmplx_t* a, const int* lda, Dblas_cmplx_t* b, const int* ldb, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRV - solve triangular linear system with single vector                       //
///////////////////////////////////////////////////////////////////////////////////////////
#define STRSV BLAS_SUBR(STRSV,strsv)
void STRSV(const char* uplo, const char* trans, const char* diag, const int* n,
const float* a, const int* lda, float* b, const int* incb, int* info);

#define DTRSV BLAS_SUBR(DTRSV,dtrsv)
void DTRSV(const char* uplo, const char* trans, const char* diag, const int* n,
const double* a, const int* lda, double* b, const int* incb, int* info);

#define CTRSV BLAS_SUBR(CTRSV,ctrsv)
void CTRSV(const char* uplo, const char* trans, const char* diag, const int* n,
const Sblas_cmplx_t* a, const int* lda, Sblas_cmplx_t* b, const int* incb, int* info);

#define ZTRSV BLAS_SUBR(ZTRSV,ztrsv)
void ZTRSV(const char* uplo, const char* trans, const char* diag, const int* n,
const Dblas_cmplx_t* a, const int* lda, Dblas_cmplx_t* b, const int* incb, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRM - solve triangular linear system with rhs matrix                          //
///////////////////////////////////////////////////////////////////////////////////////////
#define STRSM BLAS_SUBR(STRSM,strsm)
void STRSM(const char* side, const char* uplo, const char* trans, const char* diag, const int* m, const int* n,
const float* alpha, const float* a, const int* lda, float* b, const int* ldb, int* info);

#define DTRSM BLAS_SUBR(DTRSV,dtrsm)
void DTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const int* m, const int* n,
const double* alpha, const double* a, const int* lda, double* b, const int* ldb, int* info);

#define CTRSM BLAS_SUBR(CTRSV,ctrsm)
void CTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const int* m, const int* n,
const Sblas_cmplx_t* alpha, const Sblas_cmplx_t* a, const int* lda, Sblas_cmplx_t* b, const int* ldb, int* info);

#define ZTRSM BLAS_SUBR(ZTRSV,ztrsm)
void ZTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const int* m, const int* n,
const Dblas_cmplx_t* alpha, const Dblas_cmplx_t* a, const int* lda, Dblas_cmplx_t* b, const int* ldb, int* info);



///////////////////////////////////////////////////////////////////////////////////////////
//      XLARTG - compute givens rotation
///////////////////////////////////////////////////////////////////////////////////////////
#define SLARTG LAPACK_SUBR(SLARTG,slartg)
void SLARTG(const float *f, const float *g, float* cs, float* sn, float* r);

#define DLARTG LAPACK_SUBR(DLARTG,dlartg)
void DLARTG(const double *f, const double *g, double* cs, double* sn, double* r);

#define CLARTG LAPACK_SUBR(CLARTG,clartg)
void CLARTG(const Sblas_cmplx_t *f, const Sblas_cmplx_t *g, float* cs, Sblas_cmplx_t* sn, Sblas_cmplx_t* r);

#define ZLARTG LAPACK_SUBR(ZLARTG,zlartg)
void ZLARTG(const Dblas_cmplx_t *f, const Dblas_cmplx_t *g, double* cs, Dblas_cmplx_t* sn, Dblas_cmplx_t* r);

#ifdef __cplusplus
}
#endif

#endif
