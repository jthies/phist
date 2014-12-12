//! orthogonalize an mvec against an already orthogonal one. \ingroup core

//! This is the main orthogonalization routine in PHIST.           
//! It takes an orthogonal basis V and a set of vectors W,         
//! and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.         
//! (Q overwrites W here).                                         
//! The matrices R1 and R2 must be pre-allocated by the caller.    
//!                                                                
//! The algorithm used is up to numSweeps steps of classical       
//! block Gram-Schmidt, alternated with QR factorizations to       
//! normalize W. The method stops if the reduction in norm of W    
//! by a CGS step is less than approx. a factor sqrt(2).           
//!                                                                
//! If we find that W-V*R2 does not have full column rank,         
//! the matrix Q is augmented with random vectors which are made   
//! mutually orthogonal and orthogonal against V. In this case the 
//! dimension of the null space of W-V*R2 is returned in ierr>0.   
//!                                                                
//! The implementation makes use of the kernel function mvec_QR    
//! for the inner orthogonalization of W. If this function is not  
//! implemented by the kernel lib (i.e. returns -99), we use a     
//! fallback variant based on the SVQB algorithm. In this case,    
//! the relation Q*R1 = W-V*R2 does not hold, but Q is orthonormal 
//! and orthogonal against V anyway with the nullspace replaced    
//! by random vectors. To indicate the invalid R's, ierr=-7 is     
//! returned.                                                      
//!                                                                
//! If no random orthogonal vectors can be generated (after some   
//! tries) ierr=-8 is returned. This may indicate a problem with   
//! the random vector generator.                                   
//!                                                                
//! If the decrease in norm in one of the columns                  
//! in the last CGS sweep indicates that the algorithm has not     
//! yet converged, we return ierr=-9 to indicate that more steps   
//! may be advisable. This should not happen in practice if        
//! numSweeps>=2 ('twice is enough').                              
//!                                                                
void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps, int* ierr);
