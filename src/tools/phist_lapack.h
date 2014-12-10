#ifndef PHIST_LAPACK_H
#define PHIST_LAPACK_H

//TODO: cmake integration of lapacke
//      I think we should gradually move towards
//      using lapacke everywhere
#ifdef PHIST_HAVE_MKL
#include "mkl_lapack.h"
#include "mkl_lapacke.h"
#else
#include "lapacke.h"
#endif
#ifdef PHIST_SDMATS_ROW_MAJOR
#define SDMAT_FLAG LAPACK_ROW_MAJOR
#else
#define SDMAT_FLAG LAPACK_COL_MAJOR
#endif

// to allow calling fortran-style lapack interfaces, define a complex type for C
#ifdef PHIST_HAVE_MKL
typedef MKL_Complex8 Sblas_cmplx_t;
typedef MKL_Complex16 Dblas_cmplx_t;
typedef MKL_INT blas_idx_t;
#else
typedef struct Sblas_cmplx_t {
float re;
float im;
} Sblas_cmplx_t;

typedef struct Dblas_cmplx_t {
double re;
double im;
} Dblas_cmplx_t;

typedef int blas_idx_t;
#endif

// TODO - cmake/blas/lapack integration.
// this is a platform dependent macro, we should have CMake determine
// how to define the name of a fortran 77 routine
// NOTE: mkl_lapack.h defines a variety of options, so as long as it is
// used we're fine. The lower case/underscore variant here works for
// linux systems, typically.
#define LAPACK_SUBR(NAME,name) name ## _
#define BLAS_SUBR(NAME,name) name ## _

/* GEMM */
#define SGEMM BLAS_SUBR(SGEMM,sgemm)
#define DGEMM BLAS_SUBR(DGEMM,dgemm)
#define CGEMM BLAS_SUBR(CGEMM,cgemm)
#define ZGEMM BLAS_SUBR(ZGEMM,zgemm)
#define SSTEQR LAPACK_SUBR(SSTEQR,ssteqr)
#define DSTEQR LAPACK_SUBR(DSTEQR,dsteqr)
#define SGEES LAPACK_SUBR(SGEES,sgees)
#define DGEES LAPACK_SUBR(DGEES,dgees)
#define CGEES LAPACK_SUBR(CGEES,cgees)
#define ZGEES LAPACK_SUBR(ZGEES,zgees)
#define STRSEN LAPACK_SUBR(STRSEN,strsen)
#define DTRSEN LAPACK_SUBR(DTRSEN,dtrsen)
#define CTRSEN LAPACK_SUBR(CTRSEN,ctrsen)
#define ZTRSEN LAPACK_SUBR(ZTRSEN,ztrsen)
#define STREXC LAPACK_SUBR(STREXC,strexc)
#define DTREXC LAPACK_SUBR(DTREXC,dtrexc)
#define CTREXC LAPACK_SUBR(CTREXC,ctrexc)
#define ZTREXC LAPACK_SUBR(ZTREXC,ztrexc)
#define STREVC LAPACK_SUBR(STREVC,strevc)
#define DTREVC LAPACK_SUBR(DTREVC,dtrevc)
#define CTREVC LAPACK_SUBR(CTREVC,ctrevc)
#define ZTREVC LAPACK_SUBR(ZTREVC,ztrevc)
#define STRTRS LAPACK_SUBR(STRTRS,strtrs)
#define DTRTRS LAPACK_SUBR(DTRTRS,dtrtrs)
#define CTRTRS LAPACK_SUBR(CTRTRS,ctrtrs)
#define ZTRTRS LAPACK_SUBR(ZTRTRS,ztrtrs)
#define STRSV BLAS_SUBR(STRSV,strsv)
#define DTRSV BLAS_SUBR(DTRSV,dtrsv)
#define CTRSV BLAS_SUBR(CTRSV,ctrsv)
#define ZTRSV BLAS_SUBR(ZTRSV,ztrsv)
#define STRSM BLAS_SUBR(STRSM,strsm)
#define DTRSM BLAS_SUBR(DTRSV,dtrsm)
#define CTRSM BLAS_SUBR(CTRSV,ctrsm)
#define ZTRSM BLAS_SUBR(ZTRSV,ztrsm)
#define SLARTG LAPACK_SUBR(SLARTG,slartg)
#define DLARTG LAPACK_SUBR(DLARTG,dlartg)
#define CLARTG LAPACK_SUBR(CLARTG,clartg)
#define ZLARTG LAPACK_SUBR(ZLARTG,zlartg)

#ifdef PHIST_SDMATS_ROW_MAJOR
/* we might use LAPACKE for this case, but we don't really need the support
for row-major sdMats right now, I think
*/
#warning "standard lapack calls will not work for row-major sdMats"
#endif

#ifdef PHIST_HAVE_MKL
#include "mkl_blas.h"
#include "mkl_lapack.h"
#else

// we provide some C interfaces to lapack routines, it may be
// a better idea to use LAPACKE everywhere (see comment above).
// We add lapack subroutines as we go along implementing things.

#ifdef __cplusplus
extern "C" {
#endif

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEMM - general dense matrix-matrix multiplicaiton                               //
///////////////////////////////////////////////////////////////////////////////////////////
void SGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const float*, const float*, const blas_idx_t*, const float*, const blas_idx_t*, const float*,
    float *, const blas_idx_t*, int*);
void DGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const double*, const double*, const blas_idx_t*, const double*, const blas_idx_t*, const double*,
    double *, const blas_idx_t*, int*);
void CGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const Sblas_cmplx_t*, const Sblas_cmplx_t*, const blas_idx_t*, const Sblas_cmplx_t*, const blas_idx_t*, const Sblas_cmplx_t*,
    Sblas_cmplx_t *, const blas_idx_t*, int*);
void ZGEMM(const char*, const char*, const blas_idx_t*, const blas_idx_t*, const blas_idx_t*,
    const Dblas_cmplx_t*, const Dblas_cmplx_t*, const blas_idx_t*, const Dblas_cmplx_t*, const blas_idx_t*, const Dblas_cmplx_t*,
    Dblas_cmplx_t *, const blas_idx_t*, int*);

///////////////////////////////////////////////////////////////////////////////////////////
// XSTEQR - QR decomposition of symmetric tridiagonal matrices                           //
///////////////////////////////////////////////////////////////////////////////////////////


// QR decomposition of a real symmetric tridiagonal matrix
void SSTEQR(const char*, const blas_idx_t* n, float* D, float* E, float* Z, const blas_idx_t* ldz, float* work, int* info);
// QR decomposition of a real symmetric tridiagonal matrix
void DSTEQR(const char*, const blas_idx_t* n, double* D, double* E, double* Z, const blas_idx_t* ldz, double* work, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEES - Schur decomposition                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

// Schur decomposition of a real non-symmetric matrix
void SGEES(const char*, const char*, int (*select)(float*, float*), const blas_idx_t* n, 
float* a, const blas_idx_t* lda, int* sdim, float* wr, float* wi,
 float* vs, const blas_idx_t* ldvs, float* work, const blas_idx_t* lwork, blas_idx_t* bwork, int* info);
// Schur decomposition of a real non-symmetric matrix
void DGEES(const char*, const char*, int (*select)(double*, double*), const blas_idx_t* n, 
double* a, const blas_idx_t* lda, blas_idx_t*sdim, double* wr, double* wi,
 double* vs, const blas_idx_t* ldvs, double* work, const blas_idx_t* lwork, blas_idx_t* bwork, int* info);


// Schur decomposition of a complex non-hermitian matrix 
void CGEES(const char*, const char*, int (*select)(Sblas_cmplx_t*), 
const blas_idx_t* n, Sblas_cmplx_t* a, const blas_idx_t* lda, blas_idx_t* sdim,
Sblas_cmplx_t* w, Sblas_cmplx_t* vs, const blas_idx_t* ldvs, Sblas_cmplx_t* work, 
const blas_idx_t* lwork, float* rwork, blas_idx_t* bwork, int* info);
// Schur decomposition of a complex non-hermitian matrix 
void ZGEES(const char*, const char*, int (*select)(Dblas_cmplx_t*), 
const blas_idx_t* n, Dblas_cmplx_t* a, const blas_idx_t* lda, blas_idx_t* sdim,
Dblas_cmplx_t* w, Dblas_cmplx_t* vs, const blas_idx_t* ldvs, Dblas_cmplx_t* work, 
const blas_idx_t* lwork, double* rwork, blas_idx_t* bwork, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRSEN - reorder Schur form                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

//! reorder real Schur decomposition
void STRSEN(const char* job, const char* jobvs, const int *select, const int *n, float *t, 
const int *ldt, float *q, const int *ldq, float *wr, float *wi, int *m, float *s, float *sep, 
float *work, const int *lwork, int *iwork, const int *liwork, int *info);

//! reorder real Schur decomposition
void DTRSEN(const char* job, const char* jobvs, const int *select, const int *n, double *t, 
const int *ldt, double *q, const int *ldq, double *wr, double *wi, int *m, double *s, double *sep, 
double *work, const int *lwork, int *iwork, const int *liwork, int *info);

//! reorder complex Schur decomposition
void CTRSEN(const char* job, const char* jobvs, const int *select, const int *n, 
Sblas_cmplx_t *T, const int *ldt, Sblas_cmplx_t *Q, const int *ldq, Sblas_cmplx_t *W, 
const int *m, float *S, float *sep, Sblas_cmplx_t *work, const int *lwork, int *info);

//! reorder complex Schur decomposition
void ZTRSEN(const char* job, const char* jobvs, const int *select, const int *n, 
Dblas_cmplx_t *T, const int *ldt, Dblas_cmplx_t *Q, const int *ldq, Dblas_cmplx_t *W, 
const int *m, double *S, double *sep, Dblas_cmplx_t *work, const int *lwork, int *info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTREXC - reorder Schur form                                                      //
///////////////////////////////////////////////////////////////////////////////////////////

//! reorder real Schur decomposition
void STREXC(const char* compq, const int *n, float *t, const int *ldt, float *q, const int *ldq, blas_idx_t* ifst, int *ilst, float *work, int *info);

//! reorder real Schur decomposition
void DTREXC(const char* compq, const int *n, double *t, const int *ldt, double *q, const int *ldq, blas_idx_t* ifst, int *ilst, double *work, int *info);

//! reorder complex Schur decomposition
void CTREXC(const char* compq, const int *n, Sblas_cmplx_t *t, const int *ldt, Sblas_cmplx_t *q, const int *ldq, blas_idx_t* ifst, int *ilst, int *info);

//! reorder complex Schur decomposition
void ZTREXC(const char* compq, const int *n, Dblas_cmplx_t *t, const int *ldt, Dblas_cmplx_t *q, const int *ldq, blas_idx_t* ifst, int *ilst, int *info);


///////////////////////////////////////////////////////////////////////////////////////////
//      TREVC - eigenvalues and -vectors of the Schur form                               //
///////////////////////////////////////////////////////////////////////////////////////////

void STREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const float* T, const blas_idx_t* ldt, float* vl, const blas_idx_t* ldvl, float* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, float* work, int* info);

void DTREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const double* T, const blas_idx_t* ldt, double* vl, const blas_idx_t* ldvl, double* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, double* work, int* info);

  void CTREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const Sblas_cmplx_t* t, const blas_idx_t* ldt, Sblas_cmplx_t* vl, const blas_idx_t* ldvl, Sblas_cmplx_t* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, Sblas_cmplx_t* work, float* rwork, int* info);

void ZTREVC(const char* side, const char* howmny, blas_idx_t* select, const blas_idx_t* n, 
const Dblas_cmplx_t* T, const blas_idx_t* ldt, Dblas_cmplx_t* vl, const blas_idx_t* ldvl, Dblas_cmplx_t* vr, const blas_idx_t* ldvr, 
const blas_idx_t* mm, blas_idx_t* m, Dblas_cmplx_t* work, double* rwork, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRS - solve triangular linear system                                          //
///////////////////////////////////////////////////////////////////////////////////////////
void STRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const float* a, const blas_idx_t* lda, float* b, const blas_idx_t* ldb, int* info);

void DTRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const double* a, const blas_idx_t* lda, double* b, const blas_idx_t* ldb, int* info);

void CTRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const Sblas_cmplx_t* a, const blas_idx_t* lda, Sblas_cmplx_t* b, const blas_idx_t* ldb, int* info);

void ZTRTRS(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n, const blas_idx_t* nrhs, 
const Dblas_cmplx_t* a, const blas_idx_t* lda, Dblas_cmplx_t* b, const blas_idx_t* ldb, int* info);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRV - solve triangular linear system with single vector                       //
///////////////////////////////////////////////////////////////////////////////////////////
void STRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const float* a, const blas_idx_t* lda, float* b, const blas_idx_t* incb);

void DTRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const double* a, const blas_idx_t* lda, double* b, const blas_idx_t* incb);

void CTRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const Sblas_cmplx_t* a, const blas_idx_t* lda, Sblas_cmplx_t* b, const blas_idx_t* incb);

void ZTRSV(const char* uplo, const char* trans, const char* diag, const blas_idx_t* n,
const Dblas_cmplx_t* a, const blas_idx_t* lda, Dblas_cmplx_t* b, const blas_idx_t* incb);

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRM - solve triangular linear system with rhs matrix                          //
///////////////////////////////////////////////////////////////////////////////////////////
void STRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const float* alpha, const float* a, const blas_idx_t* lda, float* b, const blas_idx_t* ldb, int* info);

void DTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const double* alpha, const double* a, const blas_idx_t* lda, double* b, const blas_idx_t* ldb, int* info);

void CTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const Sblas_cmplx_t* alpha, const Sblas_cmplx_t* a, const blas_idx_t* lda, Sblas_cmplx_t* b, const blas_idx_t* ldb, int* info);

void ZTRSM(const char* side, const char* uplo, const char* trans, const char* diag, const blas_idx_t* m, const blas_idx_t* n,
const Dblas_cmplx_t* alpha, const Dblas_cmplx_t* a, const blas_idx_t* lda, Dblas_cmplx_t* b, const blas_idx_t* ldb, int* info);



///////////////////////////////////////////////////////////////////////////////////////////
//      XLARTG - compute givens rotation
///////////////////////////////////////////////////////////////////////////////////////////
void SLARTG(const float *f, const float *g, float* cs, float* sn, float* r);

void DLARTG(const double *f, const double *g, double* cs, double* sn, double* r);

void CLARTG(const Sblas_cmplx_t *f, const Sblas_cmplx_t *g, float* cs, Sblas_cmplx_t* sn, Sblas_cmplx_t* r);

void ZLARTG(const Dblas_cmplx_t *f, const Dblas_cmplx_t *g, double* cs, Dblas_cmplx_t* sn, Dblas_cmplx_t* r);

#ifdef __cplusplus
}
#endif

#endif
#endif
