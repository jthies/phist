#include "scamac_cblas.h"

extern void dcopy_(   const BLAS_INT *, const double *, const BLAS_INT *, double *, const BLAS_INT *);
extern void daxpy_(   const BLAS_INT *, const double *, const double *, const BLAS_INT *, double *, const BLAS_INT *);
extern void dscal_(   const BLAS_INT *, const double *, double *, const BLAS_INT *);

extern double ddot_(  const BLAS_INT *, const double *, const BLAS_INT *, const double *, const BLAS_INT *);
extern double dnrm2_( const BLAS_INT *, const double *, const BLAS_INT *);

extern void zcopy_(   const BLAS_INT *, const double complex *, const BLAS_INT *, double complex *, const BLAS_INT *);
extern void zaxpy_(   const BLAS_INT *, const double complex *, const double complex *, const BLAS_INT *, double complex *, const BLAS_INT *);
extern void zdscal_(   const BLAS_INT *, const double *, double complex *, const BLAS_INT *);

//extern double complex zdotc_( const BLAS_INT *, const double complex *, const BLAS_INT *, const double complex *, const BLAS_INT *);
extern void zdotc_(double complex *, const BLAS_INT *, const double complex *, const BLAS_INT *, const double complex *, const BLAS_INT *);
//extern double * zdotc_(const BLAS_INT *, const double complex *, const BLAS_INT *, const double complex *, const BLAS_INT *);
extern double dznrm2_(const BLAS_INT *, const double complex *, const BLAS_INT *);


double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY) {
  BLAS_INT my_N=N, my_incX=incX, my_incY=incY;
  double dot = ddot_(&my_N, X, &my_incX, Y, &my_incY);
  return dot;
}


double cblas_dnrm2(const int N, const double *X, const int incX) {
  BLAS_INT my_N=N, my_incX=incX;
  double nrm2 = dnrm2_(&my_N, X, &my_incX);
  return nrm2;
}

void cblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY) {
  BLAS_INT my_N=N, my_incX=incX, my_incY=incY;
  daxpy_(&my_N, &alpha, X, &my_incX, Y, &my_incY);
}

void cblas_dscal(const int N, const double alpha, double *X, const int incX) {
  BLAS_INT my_N=N, my_incX=incX;
  dscal_(&my_N, &alpha, X, &my_incX);
}

void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY) {
  BLAS_INT my_N=N, my_incX=incX, my_incY=incY;
  dcopy_(&my_N, X, &my_incX, Y, &my_incY);
}

void cblas_zdotc_sub(const int N, const double complex *X, const int incX,
                     const double complex *Y, const int incY, double complex *dotc) {
  BLAS_INT my_N=N, my_incX=incX, my_incY=incY;
  zdotc_(dotc, &my_N, X, &my_incX, Y, &my_incY);
  // *dotc = zdotc_(&my_N, X, &my_incX, Y, &my_incY);
}

double cblas_dznrm2(const int N, const double complex *X, const int incX) {
  BLAS_INT my_N=N, my_incX=incX;
  double nrm2 = dznrm2_(&my_N, X, &my_incX);
  return nrm2;
}

void cblas_zaxpy(const int N, const double complex *alpha, const double complex *X,
                 const int incX, double complex *Y, const int incY) {
  BLAS_INT my_N=N, my_incX=incX, my_incY=incY;
  zaxpy_(&my_N, alpha, X, &my_incX, Y, &my_incY);
}


void cblas_zdscal(const int N, const double alpha, double complex *X, const int incX) {
  BLAS_INT my_N=N, my_incX=incX;
  zdscal_(&my_N, &alpha, X, &my_incX);
}

void cblas_zcopy(const int N, const double complex *X, const int incX,
                 double complex *Y, const int incY) {
  BLAS_INT my_N=N, my_incX=incX, my_incY=incY;
  zcopy_(&my_N, X, &my_incX, Y, &my_incY);
}


