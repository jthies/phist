/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! wrap up an iterative solver as an operator, that
//! is op->apply(alpha,op,X,beta,Y) will actually 
//! approximate Y=(A-sigma*I)\X. alpha and beta are
//! not used, resp. must be alpha=1, beta=0. 
void SUBR(linearOp_wrap_solver)(TYPE(linearOp_ptr) Ainv_op,TYPE(const_sparseMat_ptr) A, 
        _ST_ shift, phist_ElinSolv method,int block_size, _MT_ tol,int maxIter,int* iflag);

