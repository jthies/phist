/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \brief A simple PCG implementation that works on several vectors simultaneously. \ingroup blockedPCG
//!
//! For more information on blocked solvers see phist_blockedgmres_decl.h
//! In contrast to our GMRES solver, PCG does not work with a state object but simply runs
//! until convergence (of one of the systems) or failure.
//!
//! \param [in] *nIter indicates the total max number of iterations allowed, maxIter, for any system.
//! \param [out] *nIter indicates the number of blocked iterations performed. 
//!
//! \param sol,rhs must have numSys columns 
//! \param Aop must be applicable to numSys columns
//! \param preconOp (if not NULL) must be applicable to numSys columns. 
//!
//! It is allowed that Op and preconOp act as a different linear operator on each column of the input vector, e.g. 
//! Op*X_j = (A-sigma_jB)X_j.
//!
//! \return On return, iflag will be set to: <br>
//! <0 if any error occurred, <br>
//!  0 if anyone converged and there was no error, <br>
//!  1 if the number of iterations was exceeded without any system converging.
//!
void SUBR( blockedPCG_iterate ) (TYPE(const_linearOp_ptr) Aop, 
                TYPE(const_linearOp_ptr) preconOp,
                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                int numSys, int *nIter, _MT_ const tol[], int* iflag);
