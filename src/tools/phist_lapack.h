#ifndef PHIST_LAPACK_H
#define PHIST_LAPACK_H

//TODO: cmake integration of lapacke
//      I think we should gradually move towards
//      using lapacke everywhere
#ifdef PHIST_HAVE_MKL
# include "mkl_lapack.h"
# include "mkl_lapacke.h"
typedef const char blas_char_t;
typedef MKL_Complex8 Sblas_cmplx_t;
typedef MKL_Complex16 Dblas_cmplx_t;
typedef MKL_INT blas_idx_t;
#else
# define lapack_complex_float s_complex_t
# define lapack_complex_double d_complex_t
# include "lapacke.h"
typedef lapack_complex_float Sblas_cmplx_t;
typedef lapack_complex_double Dblas_cmplx_t;
typedef int blas_idx_t;
typedef char blas_char_t;
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

#ifdef __cplusplus
extern "C" {
#endif

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
//      XLARTG - compute givens rotation
///////////////////////////////////////////////////////////////////////////////////////////
void SLARTG(const float *f, const float *g, float* cs, float* sn, float* r);

void DLARTG(const double *f, const double *g, double* cs, double* sn, double* 
r);

void CLARTG(const Sblas_cmplx_t *f, const Sblas_cmplx_t *g, float* cs, 
Sblas_cmplx_t* sn, Sblas_cmplx_t* r);

void ZLARTG(const Dblas_cmplx_t *f, const Dblas_cmplx_t *g, double* cs, 
Dblas_cmplx_t* sn, Dblas_cmplx_t* r);



#ifdef __cplusplus
} // extern "C" 
#endif

#endif
#endif
