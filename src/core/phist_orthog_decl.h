//! orthogonalize an mvec against an already orthogonal one. \ingroup core

//! This is the main orthogonalization routine in PHIST.           
//! It takes an orthogonal basis V and a set of vectors W,         
//! and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.         
//! (Q overwrites W here).                                         
//! The matrices R1 and R2 must be pre-allocated by the caller.    
//!                                                                
//! There is an option to supply an operator for B, in which case  
//! V should be a B-orthogonal basis and we B-orthogonalize W with 
//! respect to V, such that instead Q'BQ=I. This is currectly not  
//! implemented.                                                   
//!                                                                
//! The algorithm used is up to numSweeps steps of classical       
//! block Gram-Schmidt, alternated with QR factorizations to       
//! normalize W. The method stops if the reduction in norm of W    
//! by a CGS step is less than approx. a factor sqrt(2).           
//!                                                                
//! If we find that W-V*R2 does not have full column rank,         
//! the matrix Q is augmented with random vectors which are made   
//! mutually orthogonal and orthogonal against V. On exit, the     
//! variable rankWV will contain the original rank of the matrix   
//! [V,W] before the orthogonalization/randomization of W. If [V,W]
//! did not have full rank, *iflag=+1 is returned.                 
//!                                                                
//! The implementation makes use of the kernel function mvec_QR    
//! for the inner orthogonalization of W. If this function is not  
//! implemented by the kernel lib (i.e. returns -99), we use a     
//! fallback variant based on the SVQB algorithm. In this case,    
//! the relation Q*R1 = W-V*R2 does not hold, but Q is orthonormal 
//! and orthogonal against V anyway with the nullspace replaced    
//! by random vectors. To indicate the invalid R's, iflag=+2 is    
//! returned.                                                      
//!                                                                
//! If no random orthogonal vectors can be generated (after some   
//! tries) iflag=-8 is returned. This may indicate a problem with  
//! the random vector generator.                                   
//!                                                                
//! If the decrease in norm in one of the columns                  
//! in the last CGS sweep indicates that the algorithm has not     
//! yet converged, we return iflag=-9 to indicate that more steps  
//! may be advisable. This should not happen in practice if        
//! numSweeps>=2 ('twice is enough').                              
//!                                                                

void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_op_ptr) B,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps, int* rankVW, int* iflag);


//! for timing measurements provide access to a trilinos orthog manager
void SUBR(trili_orthog)(TYPE(const_mvec_ptr) V,
                        TYPE(mvec_ptr) W,
                        TYPE(const_op_ptr) B,
                        TYPE(sdMat_ptr) R1,
                        TYPE(sdMat_ptr) R2,
                        int numSweeps, int* rankVW, int* iflag);
