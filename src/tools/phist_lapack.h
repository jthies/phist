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

typedef int blas_idx_t;

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
void SGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const float*, const float*, const blas_idx_t*, const float*, const blas_idx_t*, const float*,
    float *, const blas_idx_t*, int*);
#define DGEMM BLAS_SUBR(DGEMM,dgemm)
void DGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const double*, const double*, const blas_idx_t*, const double*, const blas_idx_t*, const double*,
    double *, const blas_idx_t*, int*);
#define CGEMM BLAS_SUBR(CGEMM,cgemm)
void CGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const Sblas_cmplx_t*, const Sblas_cmplx_t*, const blas_idx_t*, const Sblas_cmplx_t*, const blas_idx_t*, const Sblas_cmplx_t*,
    Sblas_cmplx_t *, const blas_idx_t*, int*);
#define ZGEMM BLAS_SUBR(ZGEMM,zgemm)
void ZGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const Dblas_cmplx_t*, const Dblas_cmplx_t*, const blas_idx_t*, const Dblas_cmplx_t*, const blas_idx_t*, const Dblas_cmplx_t*,
    Dblas_cmplx_t *, const blas_idx_t*, int*);

///////////////////////////////////////////////////////////////////////////////////////////
// XSTEQR - QR decomposition of symmetric tridiagonal matrices                           //
///////////////////////////////////////////////////////////////////////////////////////////


// QR decomposition of a real symmetric tridiagonal matrix
#define SSTEQR LAPACK_SUBR(SSTEQR,ssteqr)
void SSTEQR(const char*, const blas_idx_t* n, float* D, float* E, float* Z, const blas_idx_t* ldz, float* work, int* info);
// QR decomposition of a real symmetric tridiagonal matrix
#define DSTEQR LAPACK_SUBR(DSTEQR,dsteqr)
void DSTEQR(const char*, const blas_idx_t* n, double* D, double* E, double* Z, const blas_idx_t* ldz, double* work, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEES - Schur decomposition                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

// Schur decomposition of a real non-symmetric matrix
#define SGEES LAPACK_SUBR(SGEES,sgees)
void SGEES(const char*, const char*, int (*select)(float*, float*), const blas_idx_t* n, 
float* a, const blas_idx_t* lda, int* sdim, float* wr, float* wi,
 float* vs, const blas_idx_t* ldvs, float* work, const blas_idx_t* lwork, blas_idx_t* bwork, int* info);
// Schur decomposition of a real non-symmetric matrix
#define DGEES LAPACK_SUBR(DGEES,dgees)
void DGEES(const char*, const char*, int (*select)(double*, double*), const blas_idx_t* n, 
double* a, const blas_idx_t* lda, blas_idx_t*sdim, double* wr, double* wi,
 double* vs, const blas_idx_t* ldvs, double* work, const blas_idx_t* lwork, blas_idx_t* bwork, int* info);


// Schur decomposition of a complex non-hermitian matrix 
#define CGEES LAPACK_SUBR(CGEES,cgees)
void CGEES(const char*, const char*, int (*select)(Sblas_cmplx_t*), 
const blas_idx_t* n, Sblas_cmplx_t* a, const blas_idx_t* lda, blas_idx_t* sdim,
Sblas_cmplx_t* w, Sblas_cmplx_t* vs, const blas_idx_t* ldvs, Sblas_cmplx_t* work, 
const blas_idx_t* lwork, float* rwork, blas_idx_t* bwork, int* info);
// Schur decomposition of a complex non-hermitian matrix 
#define ZGEES LAPACK_SUBR(ZGEES,zgees)
void ZGEES(const char*, const char*, int (*select)(Dblas_cmplx_t*), 
const blas_idx_t* n, Dblas_cmplx_t* a, const blas_idx_t* lda, blas_idx_t* sdim,
Dblas_cmplx_t* w, Dblas_cmplx_t* vs, const blas_idx_t* ldvs, Dblas_cmplx_t* work, 
const blas_idx_t* lwork, double* rwork, blas_idx_t* bwork, int* info);

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
//      XTREXC - reorder Schur form                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

//! reorder real Schur decomposition
#define STREXC LAPACK_SUBR(STREXC,strexc)
void STREXC(const char* compq, const int *n, float *t, const int *ldt, float *q, const int *ldq, blas_idx_t* ifst, int *ilst, float *work, int *info);

//! reorder real Schur decomposition
#define DTREXC LAPACK_SUBR(DTREXC,dtrexc)
void DTREXC(const char* compq, const int *n, double *t, const int *ldt, double *q, const int *ldq, blas_idx_t* ifst, int *ilst, double *work, int *info);

//! reorder complex Schur decomposition
#define CTREXC LAPACK_SUBR(CTREXC,ctrexc)
void CTREXC(const char* compq, const int *n, Sblas_cmplx_t *t, const int *ldt, Sblas_cmplx_t *q, const int *ldq, blas_idx_t* ifst, int *ilst, int *info);

//! reorder complex Schur decomposition
#define ZTREXC LAPACK_SUBR(ZTREXC,ztrexc)
void ZTREXC(const char* compq, const int *n, Dblas_cmplx_t *t, const int *ldt, Dblas_cmplx_t *q, const int *ldq, blas_idx_t* ifst, int *ilst, int *info);


///////////////////////////////////////////////////////////////////////////////////////////
//      TREVC - eigenvalues and -vectors of the Schur form                               //
///////////////////////////////////////////////////////////////////////////////////////////

#define STREVC LAPACK_SUBR(STREVC,strevc)
void STREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const float* T, const blas_idx_t* ldt, float* vl, const blas_idx_t* ldvl, float* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, float* work, int* info);

#define DTREVC LAPACK_SUBR(DTREVC,dtrevc)
void DTREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const double* T, const blas_idx_t* ldt, double* vl, const blas_idx_t* ldvl, double* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, double* work, int* info);

#define CTREVC LAPACK_SUBR(CTREVC,ctrevc)
  void CTREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const Sblas_cmplx_t* t, const blas_idx_t* ldt, Sblas_cmplx_t* vl, const blas_idx_t* ldvl, Sblas_cmplx_t* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, Sblas_cmplx_t* work, float* rwork, int* info);

#define ZTREVC LAPACK_SUBR(ZTREVC,ztrevc)
void ZTREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const Dblas_cmplx_t* T, const blas_idx_t* ldt, Dblas_cmplx_t* vl, const blas_idx_t* ldvl, Dblas_cmplx_t* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, Dblas_cmplx_t* work, double* rwork, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRS - solve triangular linear system                                          //
///////////////////////////////////////////////////////////////////////////////////////////
#define STRTRS LAPACK_SUBR(STRTRS,strtrs)
void STRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const float* a, const blas_idx_t* lda, float* b, const blas_idx_t* ldb, int* info);

#define DTRTRS LAPACK_SUBR(DTRTRS,dtrtrs)
void DTRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const double* a, const blas_idx_t* lda, double* b, const blas_idx_t* ldb, int* info);

#define CTRTRS LAPACK_SUBR(CTRTRS,ctrtrs)
void CTRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const Sblas_cmplx_t* a, const blas_idx_t* lda, Sblas_cmplx_t* b, const blas_idx_t* ldb, int* info);

#define ZTRTRS LAPACK_SUBR(ZTRTRS,ztrtrs)
void ZTRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const Dblas_cmplx_t* a, const blas_idx_t* lda, Dblas_cmplx_t* b, const blas_idx_t* ldb, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRV - solve triangular linear system with single vector                       //
///////////////////////////////////////////////////////////////////////////////////////////
#define STRSV BLAS_SUBR(STRSV,strsv)
void STRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const float* a, const blas_idx_t* lda, float* b, const blas_idx_t* incb, int* info);

#define DTRSV BLAS_SUBR(DTRSV,dtrsv)
void DTRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const double* a, const blas_idx_t* lda, double* b, const blas_idx_t* incb, int* info);

#define CTRSV BLAS_SUBR(CTRSV,ctrsv)
void CTRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const Sblas_cmplx_t* a, const blas_idx_t* lda, Sblas_cmplx_t* b, const blas_idx_t* incb, int* info);

#define ZTRSV BLAS_SUBR(ZTRSV,ztrsv)
void ZTRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const Dblas_cmplx_t* a, const blas_idx_t* lda, Dblas_cmplx_t* b, const blas_idx_t* incb, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRM - solve triangular linear system with rhs matrix                          //
///////////////////////////////////////////////////////////////////////////////////////////
#define STRSM BLAS_SUBR(STRSM,strsm)
void STRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const float* alpha, const float* a, const blas_idx_t* lda, float* b, const blas_idx_t* ldb, int* info);

#define DTRSM BLAS_SUBR(DTRSV,dtrsm)
void DTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const double* alpha, const double* a, const blas_idx_t* lda, double* b, const blas_idx_t* ldb, int* info);

#define CTRSM BLAS_SUBR(CTRSV,ctrsm)
void CTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const Sblas_cmplx_t* alpha, const Sblas_cmplx_t* a, const blas_idx_t* lda, Sblas_cmplx_t* b, const blas_idx_t* ldb, int* info);

#define ZTRSM BLAS_SUBR(ZTRSV,ztrsm)
void ZTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const Dblas_cmplx_t* alpha, const Dblas_cmplx_t* a, const blas_idx_t* lda, Dblas_cmplx_t* b, const blas_idx_t* ldb, int* info);



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
