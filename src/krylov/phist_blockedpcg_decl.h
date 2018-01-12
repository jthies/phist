/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \ingroup linear_solvers
//@{

//! \defgroup blockedCG blocked CG solver for symmetric or general linear systems
//@{

//!
//! a simple CG implementation that works on several vectors simultaneously.
//! For more information on blocked solvers see phist_blockedgmres_decl.h
//! In contrast to our GMRES solver, CG does not work with a state object but simply runs
//! until convergence (of one of the systems) or failure.
//!
//! On input, *nIter indicates the total max number of iterations allowed, maxIter, for any system.
//! On output, *nIter indicates the number of blocked iterations performed. 
//!
//! sol and rhs must both have numSys columns, and op and preconOp (if not NULL) must be applicable to numSys columns. 
//! It is allowed that Op and preconOp act as a different linear operator on each column of the input vector, e.g. 
//! Op*X_j = (A-sigma_jB)X_j.
//!
//! On return, iflag will be set to:
//!
//! <0 if any error occurred,
//!  0 if anyone converged and there was no error,
//!  1 if the number of iterations was exceeded without any system converging.
//!
void SUBR( blockedPCG_iterate ) (TYPE(const_linearOp_ptr) op, 
                TYPE(const_linearOp_ptr) preconOp,
                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                int numSys, int *nIter, _MT_ const tol[], int* iflag);

//! PCG interface for just one vector.
void SUBR( PCG_iterate ) (TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) preconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ const tol, int* iflag);
//@}
//@}
