/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

//! orthogonalize an mvec against an already orthogonal one. \ingroup core

//! This is the main orthogonalization routine in PHIST.           
//! It takes an orthogonal basis V and a set of vectors W,         
//! and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.         
//! (Q overwrites W here).                                         
//! The matrices R1 and R2 must be pre-allocated by the caller.    
//!                                                                
//! There is an option to supply an operator for B, in which case  
//! V should be a B-orthogonal basis and we B-orthogonalize W with 
//! respect to V, such that instead Q'BQ=I.                        
//!                                                                
//! The algorithm used is up to numSweeps steps of classical       
//! block Gram-Schmidt, alternated with QR factorizations to       
//! normalize W. The method stops if the reduction in norm of W    
//! by a CGS step is less than approx. a factor sqrt(2).           
//!                                                                
//! If the flag PHIST_ORTHOG_RANDOMIZE_NULLSPACE is given and we   
//! find that W-V*R2 does not have full column rank,               
//! the matrix Q is augmented with random vectors which are made   
//! mutually orthogonal and orthogonal against V. If the flag is   
//! omitted, the last (m+k)-*rankVW columns are zero. On exit, the 
//! variable rankWV will contain the original rank of the matrix   
//! [V,W] before the orthogonalization/randomization of W. If [V,W]
//! did not have full rank, *iflag=+1 is returned.                 
//!                                                                
//! If no random orthogonal vectors can be generated (after some   
//! tries) iflag=-8 is returned. This may indicate a problem with  
//! the random vector generator.                                   
//!                                                                
//! If the decrease in norm in one of the columns                  
//! in the last CGS sweep indicates that the algorithm has not     
//! yet converged, we return iflag=-9 to indicate that more steps  
//! may be advisable. Since the 'inner' and 'outer' loop (to make  
//! W'W=I and V'W=0, resp.) may deteriorate the result of one      
//! another, so that more than 2 iterations may be needed.         
//!                                                                
void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_linearOp_ptr) B,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps, int* rankVW, int* iflag);



//! This function is called by orthog and allows some additional tweaking.   
//! It assumes that B*W and WtW=W'*BW are already available. It also allows  
//! you to specify the desired orthogonalization tolerance (orthoEps) and    
//! threshold when to consider two columns as linearly dependent (rankTol).  
//! You may call this function with BW=W and B=NULL when orthogonalizing in  
//! the standard inner product.                                              
void SUBR(orthog_impl)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_linearOp_ptr) B,
                     TYPE(mvec_ptr) BW,
                     TYPE(sdMat_ptr) WtW,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps, int* rankVW, 
                     _MT_ rankTol, _MT_ orthoEps,
                     int* iflag);


