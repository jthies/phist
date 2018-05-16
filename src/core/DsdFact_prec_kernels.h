/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file DsdFact_prec_kernels.h
 * Implementation of some serial LAPaCK like functions with high precision
 * \author "Melven Roehrig-Zoellner" <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 *
 * This implementation can be used by any kernel lib that provides
 * PHIST_HIGH_PRECISION_KERNELS and implements sdMat_extract_err.
*/
#ifndef PHIST_DSDFACT_PREC_KERNELS_H
#define PHIST_DSDFACT_PREC_KERNELS_H

#include "phist_config.h"

#ifdef PHIST_SDMATS_ROWMAJOR
#error "functionality not implemented for row-major sdMats"
#endif

// threshold at which to call a matrix rank deficient
#define SINGTOL 1.0e-25

#ifdef __cplusplus
extern "C" {
#endif

//! \name some serial LAPaCK like functions with high precision
//! \addtogroup core
//!@{

//! calculates a possibly low rank approximation of a lower cholesky factor of an spd matrix
//! higher-precision + pivoting + stable low rank approximation
void phist_Dprec_cholesky(double *__restrict__ a, double *__restrict__ aC, phist_lidx n, phist_lidx lda, phist_lidx *perm, int *rank, int* iflag);

//! apply backward substitution with permuted upper triangular matrix to k vectors in col-major storage
void phist_Dprec_backwardSubst(const double *__restrict__ r, const double *__restrict__ rC, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        double *__restrict__ x, double *__restrict__ xC, phist_lidx k, phist_lidx ldx, int* iflag);

//! apply forward substitution with permuted transposed upper triangular matrix
void phist_Dprec_forwardSubst(const double *__restrict__ r, const double *__restrict__ rC, phist_lidx n, phist_lidx ldr, phist_lidx *p, int rank,
        double *__restrict__ x, double *__restrict__ xC, phist_lidx k, phist_lidx ldx, int* iflag);

//! given symmetric A=V'V, compute B s.t. Q=V*B is orthonormal. B is computed in-place,
//! if A is found to be numerically rank deficient, the last n-*rank columns of B will be zeroed out
//! s.t. Q has exactly rank *rank. In bi, biC we return the inverse of B.
void phist_Dprec_qb(double *__restrict__ a,  double *__restrict__ aC, 
                    double *__restrict__ bi, double *__restrict__ biC, 
                    phist_lidx n, phist_lidx lda, int *rank, int* iflag);

//!@}
                    
#ifdef __cplusplus
} //extern "C"
#endif

#endif
