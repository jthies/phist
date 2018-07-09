/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_linear_solvers_decl.h 
//! \brief Iterative methods for linear systems

//! \ingroup linear_solvers
//!@{

//! \brief BiCGStab for a single general non-Hermitian linear systems
void SUBR(BiCGStab)(TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) rightPreconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ tol, int* iflag);

//! \brief Conjugate Gradient method for a single Hermitian and positive definite (hpd) linear systems
void SUBR(PCG)(TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) preconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ tol, int* iflag);

//! \brief restarted GMRES implementation that may work on block_size vectors simultaneously,
//! building a separate Krylov subspace for each of them. The total number of systems to be solved
//! (num_sys) may be larger than block_size.
//!
//! High-level user interface that does not require knowledge of the state object. This interface
//! should be used if you only want run a restarted GMRES solver on several systems.
//!
//! \param [in] Aop linear operator for which to solve Aop*sol=rhs
//! \param [in] Pop  (right) preconditioner (optional)
//! \param [in] rhs right-hand side(s)
//! \param [in] sol starting guess
//! \param [out] sol solution vector(s)
//! \param [in] *nIter indicates the total max number of iterations allowed, maxIter, for any system.
//! \param [out] *nIter indicates the number of GMRES iterations.
//! \param [in] tol[] (dimension equal to number of vectors in sol and rhs), the desired tolerance for each system
//! \param max_blocks the maximum number of iterations before a restart
//!
//!
extern "C" void SUBR( restartedBlockedGMRES ) ( TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(mvec_ptr) rhs, TYPE(mvec_ptr) sol_in, int num_sys,
        int nIter[], _MT_ const tol[], int block_size, int max_blocks, int* iflag)

//! \brief restarted GMRES (GMRES(m)) for num_sys non-Hermitian linear systems.

//! This is a simplified interface to restartedBlockedGmres where num_sys==block_size and the same
//! tol is applied for each system. The number of systems num_sys is determined by the number of vectors
//! in sol and rhs.
void SUBR(restartedGMRES)( TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(mvec_ptr) rhs, TYPE(mvec_ptr) sol,
        int *nIter, _MT_ tol, int m, int* iflag);

//! \brief run the CARP-CG solver for general systems with small diagonal elements

//! CARP-CG (Gordon & Gordon 2010) is a block parallel variant of the CGMN algorithm 
//! (SOR on AA^T, accelerated by CG)
void SUBR(carp_cg)( TYPE(const_linearOp_ptr) Aop,
        TYPE(mvec_ptr) rhs, TYPE(mvec_ptr) sol,
        int *nIter, _MT_ tol, int* iflag);


//!@}
