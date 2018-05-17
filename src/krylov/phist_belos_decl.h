/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \brief The Belos package from Trilinos provides block Krylov solvers.
//! We currently provide access to four diffrent methods. \ingroup krylov
//!
//!
//! method=0: block GMRES <br>
//! method=1: pseudo-block GMRES <br>
//! method=2: (block-)conjugate gradients (CG) for SPD matrices <br>
//! method=3: pseudo-block CG
//!
//! AX=B is solved for X with any number of vectors in X and B.
//!
//!  until num_iters is reached or     
//! ||r||/||r0||<tol is achieved. *num_iters is overwritten
//! by the actual number of iterations performed.
//!
//! GMRES: at most max_blocks blocks of vectors are
//! generated. A block consists of as many vectors as there
//! are columns in X and B. If max_blocks is reached, the  
//! method is restarted. For CG, max_blocks is ignored.
//!
//! The 'pseudo' block solvers form k indep. subspaces. Equivalent
//! to running a single solver for each rhs. *nconv indicates
//! (on entry) the minimum number of systems that should have converged
//! before returning from this function. If nConv==NULL, nConv=k is assumed.
void SUBR(belos)(TYPE(const_linearOp_ptr) Op, 
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B,
        TYPE(const_linearOp_ptr) Prec, 
        _MT_ tol,int* num_iters, int max_blocks,
        int variant, int *nConv,
        int* iflag);
