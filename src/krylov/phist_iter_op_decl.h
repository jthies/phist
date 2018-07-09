/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_iter_op_decl.h 
//! \brief wrap up an iterative solver as an operator

//! \brief wrap up an iterative solver as an operator \ingroup krylov
//!
//! that is op->apply(alpha,op,X,beta,Y) will actually 
//! approximate Y=(A+shift*B)\ X. alpha and beta are
//! not used, resp. must be alpha=1, beta=0. 
//! If shift!=0 and B==NULL, B=I is assumed (the identity matrix).
//! The operator P may be NULL or represent a preconditioner for the iterative
//! scheme to be used (only used by Krylov methods such as phist_PCG, phist_GMRES etc.)
void SUBR(linearOp_wrap_solver)(TYPE(linearOp_ptr) Ainv_op,TYPE(const_sparseMat_ptr) A, 
        _ST_ shift, TYPE(const_sparseMat_ptr) B, TYPE(const_linearOp_ptr) P,
        phist_ElinSolv method, _MT_ tol,int maxIter,int* iflag);

