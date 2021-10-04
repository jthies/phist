/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file mlapack_wrapper.h
 * \brief interface to multiple-precision lapack routines used in PHIST
 * 
 * \author "Jonas Thies" <j.thies@tudelft.nl>
*/

#include "phist_config.h"

#ifndef PHIST_MLAPACK_WRAPPER_H
#define PHIST_MLAPACK_WRAPPER_H

#ifdef PHIST_HAVE_MPACK_QD

#ifdef __cplusplus
extern "C" {
#endif

//! \name interface to multiple-precision lapack routines used in PHIST
//! \addtogroup core
//!@{

//! \brief symmetric eigenvalue decomposition in simulated quad precision (using function
//! Rsyevd, cf. lapack routine dsyevd) 
void phist_Drsyev(int n, double *restrict a, double *restrict aC, int lda,
                       double *restrict w, double *restrict wC, int *iflag);

//! \brief general singular value decomposition of an m times n matrix 
void phist_Drgesvd(const char *jobu, const char *jobvt, int m, int n,
            double *restrict a, double *restrict aC, int lda, double *restrict s, double *restrict sC,
            double *restrict u, double *restrict uC, int ldu, double *restrict vt, double *restrict vtC, 
            int ldvt, int *iflag);             

//!@}
            
#ifdef __cplusplus
} //extern "C"
#endif

#endif

#endif
