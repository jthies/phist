/*! \file sdmat_syevd_prec.h
 * interface to multiple-precision lapack routine Rsyevd
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 *
*/

#include "phist_config.h"

#ifndef PHIST_SDMAT_SYEVD_PREC_H
#define PHIST_SDMAT_SYEVD_PREC_H

#ifdef PHIST_HAVE_MPACK_QD

#ifdef __cplusplus
extern "C" {
#endif

//! symmetric eigenvalue decomposition in simulated quad precision (using function
//! Rsyevd, cf. lapack routine dsyevd)
void phist_Drsyev(int n, double *restrict a, double *restrict aC, int lda,
                       double *restrict w, double *restrict wC, int *iflag);

#ifdef __cplusplus
} //extern "C"
#endif

#endif

#endif
