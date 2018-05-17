/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \brief A simple Lanczos process to compute the largest and smallest eigenvalue 
//! of a linear operator. \ingroup krylov
//!
//! This code is mainly intended as a testbed 
//! for fault tolerance (FT) capabilities in PHIST right now. Only block size 1 
//! is implemented.
//!
//! \param [in] A linearOp whose largest eigenpair is sought
//! \param [in] *numIter maximum number of iterations allowed
//!
//! \param [out] *lambda_min,*lambda_max approximations to the largest and smallest eigenvalue
//! \param [out] *numIter number of iterations performed
//!
void SUBR(simple_lanczos)(TYPE(const_linearOp_ptr) A_op,
        _MT_* lambda_min, _MT_* lambda_max, int *numIter, int* iflag);



