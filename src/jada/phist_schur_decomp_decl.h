 //!                                                                                             
 //! this function does an in-place Schur decomposition of T(1:m,1:m) into ~T and S: T = S'*~T*S.
 //! The Ritz values appear on the diagonal of ~T (for the real case there may be 2x2 blocks     
 //! for complex conjugate pairs). The nselect and nsort flags indicate in which order they      
 //! should appear:                                                                              
 //! nselect=nsort=0: unsorted                                                                   
 //! nselect>0: the first <nselect> Ritz values (if the last one is a complex conjugate pair     
 //!            <nselect+1>) appear in any order in the upper left corner of T                   
 //! 0<nsort<nselect: in addition to moving a cluster of <nselect> eigenvalues, sort the         
 //!             first <nsort> of them in the upper left corner.                                 
 //!                                                                                             
 //! Example: if the Schur form is                                                               
 //!                                                                                             
 //!     a x x x x                                                                               
 //!     0 b c x x, and the eigenvalues are |a|<|e|<|lambda([b,c;d,b])|<|f|,                     
 //!     0 d b x x                                                                               
 //!     0 0 0 e x                                                                               
 //!     0 0 0 0 f                                                                               
 //!                                                                                             
 //! we get for nselect=3, nsort=1, which=LM:                                                    
 //!                                                                                             
 //!   f x x x x                                                                                 
 //!     b c x x                                                                                 
 //!     d b x x                                                                                 
 //!         e x                                                                                 
 //!           a                                                                                 
 //!                                                                                             
 //! so we guarantee that the largest one <nsort> is in the upper left corner,                   
 //! and that the 3 largest ones appear first in any order.                                      
 //!                                                                                             
 //! The type of the array ev is in fact _CT_ (dimension m), but we passit in as void* because   
 //! we want to be able to call the function from C and C++ alike.                               
 //!                                                                                             
 void SUBR(SchurDecomp)(_ST_* T, int ldT, _ST_* S, int ldS,
         int m, int nselect, int nsort, eigSort_t which, _MT_ tol,
         void* ev, int *iflag);
         
 //! generalized Schur Decomposition, (S,T)->(~S,~T,VS,WS) such that                            
 //! (T,S) = ( VS*~S*WS^T, VS*~T*WS^T ), with ~S upper Schur and ~T upper triangular. The       
 //! generalized Eigenvalues are returned in _CT_* ev[i]. As the lapack routine (XGGES) returns 
 //! ev=alpha/beta, beta may be found to be 0. In this case, we set *iflag=1 and ev[i]=0.
 void SUBR(GenSchurDecomp)(_ST_* S, int ldS, _ST_* T, int ldT, 
                           _ST_* VS, int ldVS, _ST_* WS, int ldWS,
                           int m, int nselect, int nsort, eigSort_t which, _MT_ tol,
                           void* ev, int* iflag);

 //! reorder multiple eigenvalues in a given (partial) Schur decomposition by the smallest 
 //! residual norm of the unprojected problem must be sorted up to nselected to work correctly!
 void SUBR(ReorderPartialSchurDecomp)(_ST_* T, int ldT, _ST_* S, int ldS,
        int m, int nselected, eigSort_t which, _MT_ tol, _MT_* resNorm, void* ev, int* permutation, int *iflag);


 //! reorder multiple eigenvalues in a given (partial) generalized Schur decomposition by the smallest 
 //! residual norm of the unprojected problem must be sorted up to nselected to work correctly!
 void SUBR(ReorderPartialGenSchurDecomp)(_ST_* S, int ldS, _ST_* T, int ldT, _ST_* VS, int ldVS, _ST_* WS, int ldWS,
        int m, int nselected, eigSort_t which, _MT_ tol, _MT_* resNorm, void* ev, int* permutation, int *iflag);

