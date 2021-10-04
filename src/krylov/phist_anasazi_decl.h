/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_anasazi_decl.h 
//! \brief access to block eigensolvers of the Anasazi package from Trilinos

//! \brief The Anasazi package from Trilinos provides block eigensolvers.
//! We currently provide access to the following methods \ingroup krylov
//!
//! block Krylov Schur  (no choice given to the user)
//! 
//! \param [in] A_op  pointer to the operator A
//! \param [in] B_op  pointer to the hpd. operator B (if B==NULL, B=I is assumed)
//!
//! \param [in] variant the "variant" parameter selects the Anasazi solver, see corresponding
//! enum phist_anasaziType
//!
//! \param v0 Can be used to pass in a start-up vector(-space) (can have any number of columns). <br>
//! v0 is assumed to be orthonormal. Our current implementation of 
//! subspacejada will use only the first column of v0 and perform Arnoldi iterations
//! unless v0 has exactly minBas columns. In this case v0 is used instead of
//! running Arnoldi iterations. <br>
//! A suitable input basis can be obtained from a previous
//! run of subspacejada by pre-allocating Q with minBas vectors.
//!
//! \param which LM, SM, LR, SR, or TARGET
//! \param symmetric Symmetry properties of the matrix
//!
//! \param [out] nIter number of iterations performed
//! \param [out] iflag    return code of the solver (0 on success, negative on error, positive on warning) 
//!   

void SUBR(anasazi)(      TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) Ainv_op,
                         TYPE(const_linearOp_ptr) B_op,  int variant,
                         TYPE(const_mvec_ptr) v0,  phist_EeigSort which,
                         _MT_ tol,                 int *nEig,
                         int* nIter,               int blockDim,
                         int numBlocks,
                         int symmetric,
                         TYPE(mvec_ptr) vX,         _ST_* eigs,
                         int* iflag);
  