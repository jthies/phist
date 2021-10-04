/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! \file phist_linear_solvers_decl.h 
//! \brief Iterative methods for linear systems - simplfied high level interfaces that may take multiple right-hand sides as an mvec.

//! \ingroup linear_solvers
//!@{

//! \brief BiCGStab for a general non-Hermitian linear system
void SUBR(BiCGStab)(TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) rightPreconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ tol, int* iflag);

//! \brief Conjugate Gradient method for a Hermitian and positive definite (hpd) linear system
void SUBR(PCG)(TYPE(const_linearOp_ptr) Op,
                                TYPE(const_linearOp_ptr) preconOp,
                                TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
                                int* nIter, _MT_ tol, int* iflag);

//! \brief restarted GMRES (GMRES(m)) for a non-Hermitian linear system.

//! This is a simplified interface to restartedBlockedGmres where num_sys==block_size and the same
//! tol is applied for each system. The number of systems num_sys is determined by the number of vectors
//! in sol and rhs.
void SUBR(restartedGMRES)( TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
        int *nIter, _MT_ tol, int m, int* iflag);

//! \brief run the CARP-CG solver for general systems (may even be non-square)

//! CARP-CG (Gordon & Gordon 2010) is a block parallel variant of the CGMN algorithm 
//! (SOR on AA^T, accelerated by CG). It is particularly useful for matrices with small
//! diagonal elements.
void SUBR(carp_cg)( TYPE(const_sparseMat_ptr) A,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol,
        int *nIter, _MT_ tol, int* iflag);


//!@}
