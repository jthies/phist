/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_LAPACK_H
#define PHIST_LAPACK_H

#include "phist_config.h"

#ifndef DOXYGEN

#include "phist_typedefs.h"
//TODO: cmake integration of lapacke
//      I think we should gradually move towards
//      using lapacke everywhere
#ifdef PHIST_HAVE_MKL
# include "mkl_lapack.h"
# include "mkl_lapacke.h"
typedef const char phist_blas_char;
typedef MKL_Complex8 phist_Sblas_cmplx;
typedef MKL_Complex16 phist_Dblas_cmplx;
typedef MKL_INT phist_blas_idx;
#else
# define lapack_complex_float phist_s_complex
# define lapack_complex_double phist_d_complex
# include "lapacke.h"
typedef lapack_complex_float phist_Sblas_cmplx;
typedef lapack_complex_double phist_Dblas_cmplx;
typedef int phist_blas_idx;
typedef char phist_blas_char;
#endif

#ifdef PHIST_SDMATS_ROW_MAJOR
#define SDMAT_FLAG LAPACK_ROW_MAJOR
#else
#define SDMAT_FLAG LAPACK_COL_MAJOR
#endif

// TODO - cmake/blas/lapack integration.
// this is a platform dependent macro, we should have CMake determine
// how to define the name of a fortran 77 routine
// NOTE: mkl_lapack.h defines a variety of options, so as long as it is
// used we're fine. The lower case/underscore variant here works for
// linux systems, typically.
#ifndef PHIST_HAVE_MKL
#define LAPACK_SUBR(NAME,name) name ## _
#else
#define LAPACK_SUBR(NAME,name) name ## _
#endif

#define BLAS_SUBR(NAME,name) name ## _

/* GEMM - matrix-matrix multiplication */
#define SGEMM BLAS_SUBR(SGEMM,sgemm)
#define DGEMM BLAS_SUBR(DGEMM,dgemm)
#define CGEMM BLAS_SUBR(CGEMM,cgemm)
#define ZGEMM BLAS_SUBR(ZGEMM,zgemm)
/* GEQRT */
#define SGEQRT LAPACK_SUBR(SGEQRT,sgeqrt)
#define DGEQRT LAPACK_SUBR(DGEQRT,dgeqrt)
#define CGEQRT LAPACK_SUBR(CGEQRT,cgeqrt)
#define ZGEQRT LAPACK_SUBR(ZGEQRT,zgeqrt)
/* GEMQRT */
#define SGEMQRT LAPACK_SUBR(SGEMQRT,sgemqrt)
#define DGEMQRT LAPACK_SUBR(DGEMQRT,dgemqrt)
#define CGEMQRT LAPACK_SUBR(CGEMQRT,cgemqrt)
#define ZGEMQRT LAPACK_SUBR(ZGEMQRT,zgemqrt)
/* STEQR */
#define SSTEQR LAPACK_SUBR(SSTEQR,ssteqr)
#define DSTEQR LAPACK_SUBR(DSTEQR,dsteqr)
#define CSTEQR LAPACK_SUBR(CSTEQR,csteqr)
#define ZSTEQR LAPACK_SUBR(ZSTEQR,zsteqr)
/* GEES - compute Schur form and eigenvalues A */
#define SGEES LAPACK_SUBR(SGEES,sgees)
#define DGEES LAPACK_SUBR(DGEES,dgees)
#define CGEES LAPACK_SUBR(CGEES,cgees)
#define ZGEES LAPACK_SUBR(ZGEES,zgees)
/* GGES - compute generalized Schur form and eigenvalues of (A,B) */
#define SGGES LAPACK_SUBR(SGGES,sgges)
#define DGGES LAPACK_SUBR(DGGES,dgges)
#define CGGES LAPACK_SUBR(CGGES,cgges)
#define ZGGES LAPACK_SUBR(ZGGES,zgges)
/* TRSEN - sort Schur form */
#define STRSEN LAPACK_SUBR(STRSEN,strsen)
#define DTRSEN LAPACK_SUBR(DTRSEN,dtrsen)
#define CTRSEN LAPACK_SUBR(CTRSEN,ctrsen)
#define ZTRSEN LAPACK_SUBR(ZTRSEN,ztrsen)
/* TGSEN - sort generalized Schur form */
#define STGSEN LAPACK_SUBR(STGSEN,stgsen)
#define DTGSEN LAPACK_SUBR(DTGSEN,dtgsen)
#define CTGSEN LAPACK_SUBR(CTGSEN,ctgsen)
#define ZTGSEN LAPACK_SUBR(ZTGSEN,ztgsen)
/* TREXC - swap rows in Schur form */
#define STREXC LAPACK_SUBR(STREXC,strexc)
#define DTREXC LAPACK_SUBR(DTREXC,dtrexc)
#define CTREXC LAPACK_SUBR(CTREXC,ctrexc)
#define ZTREXC LAPACK_SUBR(ZTREXC,ztrexc)
/* TGEXC - swap rows in generalized Schur form */
#define STGEXC LAPACK_SUBR(STGEXC,stgexc)
#define DTGEXC LAPACK_SUBR(DTGEXC,dtgexc)
#define CTGEXC LAPACK_SUBR(CTGEXC,ctgexc)
#define ZTGEXC LAPACK_SUBR(ZTGEXC,ztgexc)
/* TREVC - compute eigenvectors */
#define STREVC LAPACK_SUBR(STREVC,strevc)
#define DTREVC LAPACK_SUBR(DTREVC,dtrevc)
#define CTREVC LAPACK_SUBR(CTREVC,ctrevc)
#define ZTREVC LAPACK_SUBR(ZTREVC,ztrevc)
/* TRTRS - triangular solve */
#define STRTRS LAPACK_SUBR(STRTRS,strtrs)
#define DTRTRS LAPACK_SUBR(DTRTRS,dtrtrs)
#define CTRTRS LAPACK_SUBR(CTRTRS,ctrtrs)
#define ZTRTRS LAPACK_SUBR(ZTRTRS,ztrtrs)
/* TRSV */
#define STRSV BLAS_SUBR(STRSV,strsv)
#define DTRSV BLAS_SUBR(DTRSV,dtrsv)
#define CTRSV BLAS_SUBR(CTRSV,ctrsv)
#define ZTRSV BLAS_SUBR(ZTRSV,ztrsv)
/* TRSM */
#define STRSM BLAS_SUBR(STRSM,strsm)
#define DTRSM BLAS_SUBR(DTRSM,dtrsm)
#define CTRSM BLAS_SUBR(CTRSM,ctrsm)
#define ZTRSM BLAS_SUBR(ZTRSM,ztrsm)
/* LARTG */
#define SLARTG LAPACK_SUBR(SLARTG,slartg)
#define DLARTG LAPACK_SUBR(DLARTG,dlartg)
#define CLARTG LAPACK_SUBR(CLARTG,clartg)
#define ZLARTG LAPACK_SUBR(ZLARTG,zlartg)
/* GESVD */
#define SGESVD LAPACK_SUBR(SGESVD,sgesvd)
#define DGESVD LAPACK_SUBR(DGESVD,dgesvd)
#define CGESVD LAPACK_SUBR(CGESVD,cgesvd)
#define ZGESVD LAPACK_SUBR(ZGESVD,zgesvd)

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

#ifdef __cplusplus
extern "C" {
#endif

///////////////////////////////////////////////////////////////////////////////////////////
//      XTRTRV - solve triangular linear system with single vector                       //
///////////////////////////////////////////////////////////////////////////////////////////
void STRSV(const char* uplo, const char* trans, const char* diag, const phist_blas_idx* n,
const float* a, const phist_blas_idx* lda, float* b, const phist_blas_idx* incb);

void DTRSV(const char* uplo, const char* trans, const char* diag, const phist_blas_idx* n,
const double* a, const phist_blas_idx* lda, double* b, const phist_blas_idx* incb);

void CTRSV(const char* uplo, const char* trans, const char* diag, const phist_blas_idx* n,
const phist_Sblas_cmplx* a, const phist_blas_idx* lda, phist_Sblas_cmplx* b, const phist_blas_idx* incb);

void ZTRSV(const char* uplo, const char* trans, const char* diag, const phist_blas_idx* n,
const phist_Dblas_cmplx* a, const phist_blas_idx* lda, phist_Dblas_cmplx* b, const phist_blas_idx* incb);

///////////////////////////////////////////////////////////////////////////////////////////
//      XLARTG - compute givens rotation
///////////////////////////////////////////////////////////////////////////////////////////
void SLARTG(const float *f, const float *g, float* cs, float* sn, float* r);

void DLARTG(const double *f, const double *g, double* cs, double* sn, double*
r);

void CLARTG(const phist_Sblas_cmplx *f, const phist_Sblas_cmplx *g, float* cs,
phist_Sblas_cmplx* sn, phist_Sblas_cmplx* r);

void ZLARTG(const phist_Dblas_cmplx *f, const phist_Dblas_cmplx *g, double* cs,
phist_Dblas_cmplx* sn, phist_Dblas_cmplx* r);

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEMM - matrix multiplication
///////////////////////////////////////////////////////////////////////////////////////////
void SGEMM(const char* transA, const char* transB, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
const float* alpha, const float* a,  const phist_blas_idx* lda, float* b, const phist_blas_idx* ldb, const float* beta, float* c, const phist_blas_idx* ldc);

void DGEMM(const char* transA, const char* transB, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
const double* alpha, const double* a,  const phist_blas_idx* lda, double* b, const phist_blas_idx* ldb, const double* beta, double* c, const phist_blas_idx* ldc);

void CGEMM(const char* transA, const char* transB, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
const phist_Sblas_cmplx* alpha, const phist_Sblas_cmplx* a,  const phist_blas_idx* lda, phist_Sblas_cmplx* b, const phist_blas_idx* ldb, const phist_Sblas_cmplx* beta, phist_Sblas_cmplx* c, const phist_blas_idx* ldc);

void ZGEMM(const char* transA, const char* transB, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
const phist_Dblas_cmplx* alpha, const phist_Dblas_cmplx* a,  const phist_blas_idx* lda, phist_Dblas_cmplx* b, const phist_blas_idx* ldb, const phist_Dblas_cmplx* beta, phist_Dblas_cmplx* c, const phist_blas_idx* ldc);
/*
///////////////////////////////////////////////////////////////////////////////////////////
//      XGEQRT - QR factorization
///////////////////////////////////////////////////////////////////////////////////////////

void SGEQRT( const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* nb, float* a, const phist_blas_idx* lda,
                   float* tau, const phist_blas_idx* ldT, float* work, phist_blas_idx *info);
void DGEQR( const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* nb, double* a, const phist_blas_idx* lda,
                  double* tau, const phist_blas_idx* ldT, double* work, phist_blas_idx *info);
void CGEQRT( const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* nb, phist_Sblas_cmplx* a, const phist_blas_idx* lda,
                   phist_Sblas_cmplx* tau, const phist_blas_idx* ldT, phist_Sblas_cmplx* work, phist_blas_idx *info );
void ZGEQRT( const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* nb, phist_Dblas_cmplx* a, const phist_blas_idx* lda,
                   phist_Dblas_cmplx* tau, const phist_blas_idx* ldT, phist_Dblas_cmplx* work, phist_blas_idx *info );

///////////////////////////////////////////////////////////////////////////////////////////
//      XGEMQR - 'multiply by Q' after XGEQR
///////////////////////////////////////////////////////////////////////////////////////////

void SGEMQRT( phist_blas_char* side, phist_blas_char* trans, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
             const phist_blas_idx* nb, const float* A, const phist_blas_idx* lda, const float* tau, const phist_blas_idx *ldT,
             float* c, const phist_blas_idx *ldc, float* work, phist_blas_idx* info);

void DGEMQRT( phist_blas_char* side, phist_blas_char* trans, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
             const phist_blas_lidx* nb, const double* A, const phist_blas_idx* lda, const double* tau, const phist_blas_idx *ldT,
             double* c, const phist_blas_idx *ldc, double* work, phist_blas_idx* info);

void CGEMQRT( phist_blas_char* side, phist_blas_char* trans, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
             const phist_blas_lidx* nb, const phist_Sblas_cmplx* A, const phist_blas_idx* lda, const phist_Sblas_cmplx* tau, const phist_blas_idx *ldT,
             phist_Sblas_cmplx* c, const phist_blas_idx *ldc, phist_Sblas_cmplx* work, phist_blas_idx* info);

void ZGEMQRT( phist_blas_char* side, phist_blas_char* trans, const phist_blas_idx* m, const phist_blas_idx* n, const phist_blas_idx* k,
             const phist_blas_idx* nb, const phist_Dblas_cmplx* A, const phist_blas_idx* lda, const phist_Dblas_cmplx* tau, const phist_blas_idx *ldT, 
             phist_Dblas_cmplx* c, const phist_blas_idx *ldc, phist_Dblas_cmplx* work, phist_blas_idx* info);

*/
#ifdef __cplusplus
} // extern "C"
#endif

#endif

#endif /* DOXYGEN */
#endif /* PHIST_LAPACK_H */
