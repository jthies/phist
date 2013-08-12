//                                                               
// This is the main orthogonalization routine in PHIST.          
// It takes an orthogonal basis V and a set of vectors W,        
// and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.        
// (Q overwrites W here).                                        
// The matrices R1 and R2 must be pre-allocated by the caller.   
//                                                               
// The algorithm used is up to numSweeps steps of classical      
// block Gram-Schmidt, alternated with QR factorizations to      
// normalize W. The method stops if the reduction in norm of W   
// by a CGS step is less than approx. a factor sqrt(2).          
//                                                               
// If we find that W does not have full column rank,             
// the matrix Q is augmented with random vectors which are made  
// mutually orthogonal and orthogonal against V. The original    
// rank of W is returned in ierr at the end of the routine. If   
// it happens somewhere during the process, we return an error   
// code ierr=-7.                                                 
//                                                               
// If a breakdown occurs, indicating that one of the columns of  
// W lives in the space spanned by the columns of V, ierr=-8 is  
// returned. A more convenient behavior may be added later, like 
// randomizing the column(s) as before.                          
//                                                               
// If the decrease in norm in one of the columns                 
// in the last CGS sweep indicates that the algorithm has not    
// yet converged, we return ierr=-9 to indicate that more steps  
// may be advisable. This should not happen in practice if       
// numSweeps>=2 ('twice is enough').                             
//                                                               
void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps, int* ierr);
