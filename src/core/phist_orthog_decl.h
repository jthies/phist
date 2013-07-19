//! This is the main orthogonalization routine in PHIST.          
//! It takes an orthogonal basis V and a set of vectors W,        
//! and computes [Q,R1,R2] such that Q*R1 = W-V*R2, Q'Q=I.        
//! (Q overwrites W here).                                        
//! The matrices R1 and R2 must be pre-allocated by the caller.   
//! If relax=0, we guarantee orthogonality to machine precision.  
//! If not, we just check and report an orthogonalization failure 
//! as ierr=1.                                                    
void _SUBR_(orthog)(_TYPE_(const_mvec_ptr) V,
                     _TYPE_(mvec_ptr) W,
                     _TYPE_(sdMat_ptr) R1,
                     _TYPE_(sdMat_ptr) R2,
                     int relax, int* ierr); 
