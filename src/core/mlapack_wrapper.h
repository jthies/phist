/*! \file mlapack_wrapper.h
 * interface to multiple-precision lapack routines used in PHIST
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"

#ifndef PHIST_MLAPACK_WRAPPER_H
#define PHIST_MLAPACK_WRAPPER_H

#ifdef PHIST_HAVE_MPACK_QD

#ifdef __cplusplus
extern "C" {
#endif

//! symmetric eigenvalue decomposition in simulated quad precision (using function
//! Rsyevd, cf. lapack routine dsyevd)
void phist_Drsyev(int n, double *restrict a, double *restrict aC, int lda,
                       double *restrict w, double *restrict wC, int *iflag);

//! general singular value decomposition of an m times n matrix
void phist_Drsyev(const char *jobu, const char *jobvt, int m, int n,
            double *restrict a, double *restrict aC, int lda, double *restrict s, double *restrict sC,
            double *restrict u, double *restrict uC, int ldu, double *restrict vt, double *restrict vtC, 
            int ldvt, int *iflag);             

#ifdef __cplusplus
} //extern "C"
#endif

#endif

#endif
