/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_jadaOp_decl.h 
//! \brief routines to create operator for Jacobi-Davidson and projection operators

//! \addtogroup jada
//!@{

//! \brief Create projected and shifted operator for Jacobi-Davidson,
//! op*X = (I-BV*V')(A*X+B*X*sigma)(I-V*BV')
//! by setting the appropriate pointers.
//!
//! No data is copied.
//! For B_op==NULL, the pre-projection (I-V*BV') is
//! omitted to save reductions (and BV is ignored).
//! 
//! We make use of the apply_shifted function in the given A_op object. To implement
//! the operation (A+sigma*B)X for B!=I, the A_op should implement this in apply_shifted.
//! To construct such an operator from two sparseMats, linearOp_wrap_sparseMat_pair
//! can be used.
void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(linearOp_ptr)          jdOp,    int*                  iflag);

//! Deletes the operator jadaOp
void SUBR(jadaOp_delete)(TYPE(linearOp_ptr)  jdOp, int *iflag);

//! \brief Create a preconditioner for the inner solve in Jacobi-Davidson.                                      
//!                                                                                                      
//! Given a linear operator that is a preconditioner for A-sigma_j*B, this function will simply          
//! wrap it up to use apply_shifted when apply() is called. We need this because our implementations     
//! of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-    
//! ditioning operator. Obviously not all preconditioners are able to handle varying shifts without      
//! recomputing, this is not taken into account by this function:in that case the input P_op must be     
//! updated beforehand, or the existing preconditioner for e.g. A or A-tau*B is applied.                 
//!                                                                                                      
//! If V is given, the preconditioner application will include either <br>
//!
//! a skew-projection (if projType==1) <br>
//! Y <- (I - P_op\ V (BV'P_op\ V)^{-1} (BV)') P_op\ X                                                      
//!                                                                                                      
//! a regular projection (if projType==0) <br>
//! Y <- (I - V(BV)') P_op\ X                                                                             
//!                                                                                                      
//! If BV==NULL, BV=V is assumed.                                                                        
void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op, 
                           TYPE(const_mvec_ptr)  V,
                           TYPE(const_mvec_ptr)  BV,
                           const _ST_ sigma[], int nvec, 
                           TYPE(linearOp_ptr) jdPrec,
                           int projType,
                           int* iflag);

//! Deletes the preconditioner
void SUBR(jadaPrec_delete)(TYPE(linearOp_ptr) jdPrec, int *iflag);

//! \brief Add a left preconditioner created by jadaPrec_create to a jadaOp.

//! The effect of the apply function will afterwards by Y <- alpha*(jadaPrec*jadaOp*X) + beta*Y,
//! the projections used are determined by the AB_op and jadaPrec operators. If jadaPrec==NULL, 
//! the operator is reset to it's original effect.
void SUBR(jadaOp_set_leftPrecond)(TYPE(linearOp_ptr) jadaOp, TYPE(const_linearOp_ptr) jadaPrec, int* iflag); 

//! \brief Create projection operator
//!
//! Y <- alpha*(I - W*V')X + beta*Y
void SUBR(projOp_create)(TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) W, TYPE(linearOp_ptr) proj_Op, int* iflag);

//! Delete projection operator
void SUBR(projOp_delete)(TYPE(linearOp_ptr) proj_Op, int *iflag);

//! \brief Create skew projection operator for a preconditioner P_op
//!
//! Y <- (I - P_op\ V (BV'P_op\ V)^{-1} (BV)') * X
void SUBR(skewProjOp_create)(TYPE(const_linearOp_ptr) P_op,
        TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) BV,
        TYPE(linearOp_ptr) skew_Op, int* iflag);

//! \brief Create projected and shifted operator for Jacobi-Davidson using the linearOp_product operator
//! and a variable combination of projections.

//! \param method With method we can choose, which Projections to use: <br>
//! "NONE": op*X = (A+B*sigma)*X (only shifted operator), <br>
//! "PRE": op*X = (A+B*sigma)(I-V*BV')*X (with pre-projection), <br>
//! "POST": op*X = (I-BV*V')(A+B*sigma)*X (with post-projection), <br>
//! "PRE_POST": op*X = (I-BV*V')(A+B*sigma)(I-V*BV')*X (with pre- and post-projection), <br>
//! "SKEW": op*X = (I-(K\ V)*((BV)'K\ V)^{-1}*(BV)')*(K\+sigma*I)*(A+B*sigma)*X (with skew-projection and
//! preconditioner K), <br>
//! "ALL": op*X = P_{skew}*(K\+sigma*I)*P_{proj}^T*(A+B*sigma)*P_{proj}*X (with pre-, post-,
//! skew-projection and preconditioner)
//!
//! \param onlyPrec will change the methods "SKEW" and "ALL": <br>
//! onlyPrec == 0: we use skew-projection and preconditioner, <br>
//! onlyPrec == 1: we only use the preconditioner
void SUBR(jadaOp_create_impl)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op, TYPE(const_linearOp_ptr)    Prec_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], const _ST_ sigma_prec[],int     nvec,
                         TYPE(linearOp_ptr)          jdOp, phist_Eprojection method,
                         int onlyPrec, int*                  iflag);

//!@}
