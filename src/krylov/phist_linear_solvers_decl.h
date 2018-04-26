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

//! \brief BiCGStab for a single general non-Hermitian linear systems
void SUBR(BiCGStab)(TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) rightPreconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ const tol, int* iflag);

//! Conjugate Gradient method for a single Hermitian and positive definite (hpd) linear systems
void SUBR(PCG)(TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) preconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ const tol, int* iflag);

//! \brief restarted GMRES implementation that may work on several vectors simultaneously,
//! building a separate Krylov subspace for each of them.
//!
//! high-level user interface that does not require knowledge of the state object. This interface
//! should be used if you only want run a restarted GMRES solver on one or more systems
//!
//! \param [in] *nIter indicates the total max number of iterations allowed, maxIter, for any system.
//! \param [out] *nIter indicates the number of GMRES iterations.
//!
//! \param [out] sol_in gives the solution of the restarted GMRES
//!
void SUBR(restartedBlockedGMRES)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, int numSys,
        int* nIter, _MT_ const tol[], int block_size, int max_blocks, int* iflag);

//! restarted GMRES (GMRES(m)) for a single non-Hermitian linear system
void SUBR(restarteedGMRES)( TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(mvec_ptr) rhs, TYPE(mvec_ptr) sol,
        int *nIter, _MT_ tol, int m, int* iflag);


//@}
