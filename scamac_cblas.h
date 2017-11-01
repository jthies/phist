/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_CBLAS_H
#define SCAMAC_CBLAS_H

#include <complex.h>

/* basic wrapper for the few BLAS routines needed here */

/* COPYRIGHT: CBLAS .... */

#define BLAS_INT int

double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
double cblas_dnrm2(const int N, const double *X, const int incX);

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY);
void cblas_dscal(const int N, const double alpha, double *X, const int incX);
void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);

void cblas_zdotc_sub(const int N, const double complex *X, const int incX,
                     const double complex *Y, const int incY, double complex *dotc);
double cblas_dznrm2(const int N, const double complex *X, const int incX);

void cblas_zaxpy(const int N, const double complex *alpha, const double complex *X,
                 const int incX, double complex *Y, const int incY);
void cblas_zdscal(const int N, const double alpha, double complex *X, const int incX);
void cblas_zcopy(const int N, const double complex *X, const int incX,
                 double complex *Y, const int incY);


#endif /* SCAMAC_CBLAS_H */
