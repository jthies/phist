/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
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
# include "mkl_blas.h"
# include "mkl_lapacke.h"
typedef const char phist_blas_char;
typedef MKL_Complex8 phist_Sblas_cmplx;
typedef MKL_Complex16 phist_Dblas_cmplx;
typedef MKL_INT phist_blas_idx;
#define BLAS_SUBR(NAME,name) name
#define LAPACK_SUBR(NAME,name) name
#define LAPACKE_SUBR(NAME,name) LAPACKE_##name
#else
/* this works for OpenBLAS, not sure about other lapack vendors... */
# include "lapacke.h"
# ifndef lapack_complex_float
#  define lapack_complex_float phist_s_complex
#  define lapack_complex_double phist_d_complex
# endif

typedef lapack_complex_float phist_Sblas_cmplx;
typedef lapack_complex_double phist_Dblas_cmplx;
typedef int phist_blas_idx;
typedef char phist_blas_char;

#define BLAS_SUBR(NAME,name) name##_
#define LAPACK_SUBR(NAME,name) LAPACK_ ## name
#define LAPACKE_SUBR(NAME,name) LAPACKE_ ## name
#endif

#ifdef PHIST_SDMATS_ROW_MAJOR
#define SDMAT_FLAG LAPACK_ROW_MAJOR
#else
#define SDMAT_FLAG LAPACK_COL_MAJOR
#endif

/* GEMM - matrix-matrix multiplication */
#define SGEMM BLAS_SUBR(SGEMM,sgemm)
#define DGEMM BLAS_SUBR(DGEMM,dgemm)
#define CGEMM BLAS_SUBR(CGEMM,cgemm)
#define ZGEMM BLAS_SUBR(ZGEMM,zgemm)
/* GEQRT */
#define SGEQRT LAPACKE_SUBR(SGEQRT,sgeqrt)
#define DGEQRT LAPACKE_SUBR(DGEQRT,dgeqrt)
#define CGEQRT LAPACKE_SUBR(CGEQRT,cgeqrt)
#define ZGEQRT LAPACKE_SUBR(ZGEQRT,zgeqrt)
/* GEMQRT */
#define SGEMQRT LAPACKE_SUBR(SGEMQRT,sgemqrt)
#define DGEMQRT LAPACKE_SUBR(DGEMQRT,dgemqrt)
#define CGEMQRT LAPACKE_SUBR(CGEMQRT,cgemqrt)
#define ZGEMQRT LAPACKE_SUBR(ZGEMQRT,zgemqrt)
/* STEQR */
#define SSTEQR LAPACKE_SUBR(SSTEQR,ssteqr)
#define DSTEQR LAPACKE_SUBR(DSTEQR,dsteqr)
#define CSTEQR LAPACKE_SUBR(CSTEQR,csteqr)
#define ZSTEQR LAPACKE_SUBR(ZSTEQR,zsteqr)
/* GEES - compute Schur form and eigenvalues A */
#define SGEES LAPACKE_SUBR(SGEES,sgees)
#define DGEES LAPACKE_SUBR(DGEES,dgees)
#define CGEES LAPACKE_SUBR(CGEES,cgees)
#define ZGEES LAPACKE_SUBR(ZGEES,zgees)
/* GGES - compute generalized Schur form and eigenvalues of (A,B) */
#define SGGES LAPACKE_SUBR(SGGES,sgges)
#define DGGES LAPACKE_SUBR(DGGES,dgges)
#define CGGES LAPACKE_SUBR(CGGES,cgges)
#define ZGGES LAPACKE_SUBR(ZGGES,zgges)
/* TRSEN - sort Schur form */
#define STRSEN LAPACKE_SUBR(STRSEN,strsen)
#define DTRSEN LAPACKE_SUBR(DTRSEN,dtrsen)
#define CTRSEN LAPACKE_SUBR(CTRSEN,ctrsen)
#define ZTRSEN LAPACKE_SUBR(ZTRSEN,ztrsen)
/* TGSEN - sort generalized Schur form */
#define STGSEN LAPACKE_SUBR(STGSEN,stgsen)
#define DTGSEN LAPACKE_SUBR(DTGSEN,dtgsen)
#define CTGSEN LAPACKE_SUBR(CTGSEN,ctgsen)
#define ZTGSEN LAPACKE_SUBR(ZTGSEN,ztgsen)
/* TREXC - swap rows in Schur form */
#define STREXC LAPACKE_SUBR(STREXC,strexc)
#define DTREXC LAPACKE_SUBR(DTREXC,dtrexc)
#define CTREXC LAPACKE_SUBR(CTREXC,ctrexc)
#define ZTREXC LAPACKE_SUBR(ZTREXC,ztrexc)
/* TGEXC - swap rows in generalized Schur form */
#define STGEXC LAPACKE_SUBR(STGEXC,stgexc)
#define DTGEXC LAPACKE_SUBR(DTGEXC,dtgexc)
#define CTGEXC LAPACKE_SUBR(CTGEXC,ctgexc)
#define ZTGEXC LAPACKE_SUBR(ZTGEXC,ztgexc)
/* TREVC - compute eigenvectors */
#define STREVC LAPACKE_SUBR(STREVC,strevc)
#define DTREVC LAPACKE_SUBR(DTREVC,dtrevc)
#define CTREVC LAPACKE_SUBR(CTREVC,ctrevc)
#define ZTREVC LAPACKE_SUBR(ZTREVC,ztrevc)
/* TRTRS - triangular solve */
#define STRTRS LAPACKE_SUBR(STRTRS,strtrs)
#define DTRTRS LAPACKE_SUBR(DTRTRS,dtrtrs)
#define CTRTRS LAPACKE_SUBR(CTRTRS,ctrtrs)
#define ZTRTRS LAPACKE_SUBR(ZTRTRS,ztrtrs)
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
/* LARTG: C interface missing in MKL and OpenBLAS, so we switch to [S|D]LARTGP */
/* and call the Fortran variant [c|z]lartgp_ directly because that one is mis- */
/* sing in OpenBLAS, too.                                                      */
#define SLARTGP LAPACKE_SUBR(SLARTGP,slartgp)
#define DLARTGP LAPACKE_SUBR(DLARTGP,dlartgp)
#ifdef PHIST_HAVE_MKL
#define CLARTG LAPACK_SUBR(CLARTG,clartg)
#define ZLARTG LAPACK_SUBR(ZLARTG,zlartg)
#else
#define CLARTG LAPACK_GLOBAL(clartg, CLARTG)
#define ZLARTG LAPACK_GLOBAL(zlartg, ZLARTG)
# ifdef __cplusplus
extern "C" {
# endif
/* missing in OpenBLAS lapack.h */
void CLARTG(lapack_complex_float* F, lapack_complex_float* G, float* CS, lapack_complex_float* SN, lapack_complex_float* R);
void ZLARTG(lapack_complex_double* F, lapack_complex_double* G, double* CS, lapack_complex_double* SN, lapack_complex_double* R);
# ifdef __cplusplus
}
# endif
#endif
/* GESVD */
#define SGESVD LAPACKE_SUBR(SGESVD,sgesvd)
#define DGESVD LAPACKE_SUBR(DGESVD,dgesvd)
#define CGESVD LAPACKE_SUBR(CGESVD,cgesvd)
#define ZGESVD LAPACKE_SUBR(ZGESVD,zgesvd)

#ifdef PHIST_SDMATS_ROW_MAJOR
/* we might use LAPACKE for this case, but we don't really need the support
for row-major sdMats right now, I think
*/
#warning "standard lapack calls will not work for row-major sdMats"
#endif

#ifndef PHIST_HAVE_MKL

/* OpenBLAS has some conflicting declarations of BLAS routines
   in the headers f77blas.h and lapack.h (!?), so we can't include
   both. The workaround here is to declare all BLAS routines we need
   ourselves. An alternative would be to use cblas.h, but that's a
   slightly different interface, and not quite complete either in OpenBLAS.
 */

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

#ifdef __cplusplus
} // extern "C"
#endif

#endif

#endif /* DOXYGEN */
#endif /* PHIST_LAPACK_H */
