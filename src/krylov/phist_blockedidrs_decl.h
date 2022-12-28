/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_blockedbicgstab_decl.h
//! \brief blocked BiCGStab solver for symmetric or general linear systems

//! \brief IDR(s) linear solver that works on several vectors simultaneously.
//! In contrast to other 'blocked' iterative solvers, this implementation builds
//! a single Krylov space for the block diagonal matrix kron(I, Aop[i]), where 
//! Aop[i] is the effect that Aop has on column i. Convergence is measured based
//! on the complete system and achieved if ||Aop x[:] - b[:] ||_F<tol.
//! \ingroup blockedIDRS
//!
//! \param [in] *nIter indicates the total max number of iterations allowed, maxIter, for any system.
//! \param [out] *nIter indicates the number of blocked iterations performed. 
//!
//! \param sol,rhs must have numSys columns
//! \param Op must be applicable to numSys columns. It is allowed that
//! Op acts as a different linear operator on each column of the input vector, e.g. Op*X_j = (A-sigma_jB)X_j
//!
//! \param V is ignored right now.
//!
//! \param s>=1 indicates the length of the recurrence (the number of vectors against which to  explicitly orthogonalize)
//!
//! \return On return, iflag will be set to: <br>
//! <0 if any error occurred, <br>
//!  0 if anyone converged and there was no error, <br>
//!  1 if the number of iterations was exceeded without any system converging.
//!
void SUBR( blockedIDRs_iterate ) (TYPE(const_linearOp_ptr) Op,
                TYPE(const_linearOp_ptr) rightPreconOp,
                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol, TYPE(const_mvec_ptr) V,
                int numSys, int *nIter, _MT_ const tol, int s, int* iflag);
