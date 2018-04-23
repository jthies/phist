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

//! \defgroup blockedQMR blocked QMR solver for symmetric or general linear systems
//@{


//! \brief A simple QMR implementation that works on several vectors simultaneously,
//! building a separate Krylov subspace for each of them.
//!
//! For more information on blocked solvers see phist_blockedgmres_decl.h
//! In contrast to our GMRES solver, QMR does not work with a state object but simply runs
//! until convergence (of one of the systems) or failure.
//!
//! \param *nIter On input, *nIter indicates the total max number of iterations allowed, maxIter, for any system.
//! On output, *nIter indicates the number of blocked iterations performed. 
//!
//! \param sol must have numSys columns
//! \param rhs must have numSys columns
//! \param Op must be applicable to numSys columns. It is allowed that
//! Op acts as a different linear operator on each column of the input vector, e.g. Op*X_j = (A-sigma_jB)X_j
//!
//! \param V If V!=NULL is given, it represents the approximate eigenspace Q~ corresponding to the shifts inside Jacobi-Davidson
//! and can be used for checking convergence based on the eigenvalue residuals rather than the linear system residuals.
//!
//! \return On return, iflag will be set to:
//! <0 if any error occurred,
//!  0 if anyone converged and there was no error,
//!  1 if the number of iterations was exceeded without any system converging.
//!
void SUBR( blockedQMR_iterate ) (TYPE(const_linearOp_ptr) Op, 
                TYPE(const_linearOp_ptr) rightPreconOp,
                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol, TYPE(const_mvec_ptr) V,
                int numSys, int *nIter, _MT_ const tol[], int symmetry, int* iflag);

//@}
//@}
