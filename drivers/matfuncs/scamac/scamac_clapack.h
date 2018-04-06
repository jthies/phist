/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_CLAPACK_H
#define SCAMAC_CLAPACK_H

#include <complex.h>

/* basic wrapper for the few BLAS routines needed here */

/* COPYRIGHT: CBLAS .... */

extern void dsteqr_(char *compz, int *n, double *d, double *e, double *z, int *ldz, double *work, int *info);
extern void dlarnv_(int *idist, int *iseed, int *n, double *x);
extern void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork, int *info);

typedef struct {
  double re, im;
} LAPACK_COMPLEX;
extern void zheev_(char *jobz, char *uplo, int *n, LAPACK_COMPLEX *a, int *lda, double *w, LAPACK_COMPLEX *work, int *lwork, double *rwork, int *info);

void dgees_(char *jobvs, char *sort, void *select(), int *n, double *a, int *lda, int *sdim, double *wr,
            double *wi, double *vs, int *ldvs, double *work, int *lwork, int *bwork, int *info);

void zgees_(char *jobvs, char *sort, void *select(), int *n, LAPACK_COMPLEX *a, int *lda, int *sdim,
            LAPACK_COMPLEX *w, LAPACK_COMPLEX *vs, int *ldvs, LAPACK_COMPLEX *work, int *lwork,
            double *rwork, int *bwork, int *info);

#endif /* SCAMAC_CLAPACK_H */
