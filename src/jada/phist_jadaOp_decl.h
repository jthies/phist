/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! create projected and shifted operator for Jacobi-Davidson,
//! op*X = (I-BV*V')(A*X+B*X*sigma)(I-V*BV')
//! by setting the appropriate pointers. No data is copied.
//! For B_op==NULL, the pre-projection I-V*BV') is
//! omitted to save reductions (and BV is ignored).
//! 
//! We make use of the apply_shifted function in the given A_op object. To implement
//! the operation (A+sigma*B)X for B!=I, the A_op should implement this in apply_shifted.
//! To construct such an operator from two sparseMats, SUBR(linearOp_wrap_sparseMat_pair)
//! can be used.
void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(linearOp_ptr)          jdOp,    int*                  iflag);

void SUBR(jadaOp_delete)(TYPE(linearOp_ptr)  jdOp, int *iflag);

//! create a preconditioner for the inner solve in Jacobi-Davidson.                                      
//!                                                                                                      
//! Given a linear operator that is a preconditioner for A-sigma_j*B, this function will simply          
//! wrap it up to use apply_shifted when apply() is called. We need this because our implementations     
//! of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-    
//! ditioning operator. Obviously not all preconditioners are able to handle varying shifts without      
//! recomputing, this is not taken into account by this function:in that case the input P_op must be     
//! updated beforehand, or the existing preconditioner for e.g. A or A-tau*B is applied.                 
//!                                                                                                      
//! If V is given, the preconditioner application will include either                                    
//!                                                                                                      
//! a skew-projection (if projType==1)                                                                   
//!                                                                                                      
//! Y <- (I - P_op\V (BV'P_op\V)^{-1} (BV)') P_op\X                                                      
//!                                                                                                      
//! a regular projection (if projType==0)                                                                
//!                                                                                                      
//! Y <- (I - V(BV)') P_op\X                                                                             
//!                                                                                                      
//! If BV==NULL, BV=V is assumed.                                                                        
void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op, 
                           TYPE(const_mvec_ptr)  V,
                           TYPE(const_mvec_ptr)  BV,
                           const _ST_ sigma[], int nvec, 
                           TYPE(linearOp_ptr) jdPrec,
                           int projType,
                           int* iflag);

//! add a left preconditioner created by jadaPrec_create to a jadaOp.

//! The effect of the apply function will afterwards by Y <- alpha*(jadaPrec*jadaOp*X) + beta*Y,
//! the projections used are determined by the AB_op and jadaPrec operators. If jadaPrec==NULL, 
//! the operator is reset to it's original effect.
void SUBR(jadaOp_set_leftPrecond)(TYPE(linearOp_ptr) jadaOp, TYPE(const_linearOp_ptr) jadaPrec, int* iflag);

//! create pre projection operator for Jacobi-Davidson,
//! op*X = (I-V*BV')X
//! For B_op==NULL we get op*X = (I-VV')X
//! 
//! We make use of the apply function in the given B_op object.
void SUBR(pre_projection_Op_create)(TYPE(const_linearOp_ptr) B_op, 
                                    TYPE(const_mvec_ptr) V,
	                                TYPE(const_mvec_ptr) BV, 
									TYPE(linearOp_ptr) pre_proj_Op, 
									int* iflag);
									
//! create post projection operator for Jacobi-Davidson,
//! op*X = (I-BV*V')X
//! For B_op==NULL we get op*X = (I-VV')X
//! 
//! We make use of the apply function in the given B_op object.
void SUBR(post_projection_Op_create)(TYPE(const_linearOp_ptr) B_op, 
                                    TYPE(const_mvec_ptr) V,
	                                TYPE(const_mvec_ptr) BV, 
									TYPE(linearOp_ptr) pre_proj_Op, 
									int* iflag);	

// shifted operator
//for B==NULL: y <- alpha*(A + sigma*I)*X + beta*Y
//for B!=NULL: y <- alpha*(A + sigma*B)*X + beta*Y
extern "C" void SUBR(shifted_Op_create)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(linearOp_ptr)          shift_Op,    int*                  iflag);									