/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! The Anasazi package from Trilinos provides block eigensolvers.
//! We currently provide access to the following methods
//!
//! block Krylov Schur  (no choice given to the user)
//! 
//! input parameters: similar to subspacejada (TODO: documentation)
//!
//! the "variant" parameter selects the Anasazi solver, see corresponding
//! enum phist_anasaziType
void SUBR(anasazi)(      TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) Ainv_op,
                         TYPE(const_linearOp_ptr) B_op,  int variant,
                         TYPE(const_mvec_ptr) v0,  phist_EeigSort which,
                         _MT_ tol,                 int *nEig,
                         int* nIter,               int blockDim,
                         int numBlocks,
                         int symmetric,
                         TYPE(mvec_ptr) vX,         _ST_* eigs,
                         int* iflag);

