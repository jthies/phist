/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_simple_arnoldi_decl.h 
//! \brief Implementation of simple Arnoldi routines

//! \brief A simple Arnoldi process to start up the JaDa iteration. \ingroup krylov
//!
//!
//! Given a minimum basis size m, compute V(:,1:m+1), 
//! H(1:m+1,1:m) such that A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m)
//!
//! \param [in] v0,V,H allocated with m+1 resp. m columns
//! and nloc resp. m+1 rows.
//!
//! If B_op!=NULL, B is assumed to be symmetric and positive definite 
//! and the B-inner product is used so that 
//! A*V(:,1:m) = V(:,1:m+1)*H(1:m+1,1:m), and V^T*B*V=I.
//!
void SUBR(simple_arnoldi)(TYPE(const_linearOp_ptr) A_op, TYPE(const_linearOp_ptr) B_op, TYPE(const_mvec_ptr) v0,
        TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, int m, int* iflag);

//! starts with random vector! \ingroup krylov
void SUBR(simple_blockArnoldi)(TYPE(const_linearOp_ptr) A_op, TYPE(const_linearOp_ptr) B_op,
                               TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV,
                               TYPE(sdMat_ptr) H, int m, int blockSize, int* iflag);
