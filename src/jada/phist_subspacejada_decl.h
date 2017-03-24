/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \addtogroup jada
//@{

//! \defgroup bjdqr Jacobi-Davidson QR method with blocking, locking and restart
//@{


//! Subspace Jacobi-Davidson for extreme eigenvalues

//! Tries to compute a partial Schur form $(Q,R)$ of dimension opts.nEigs                       
//! of the stencil $A*x-\lambda*B*x$ with a general linear operator $A$ and a                   
//! hermitian positive definite (hpd.) linear operator $B$ using a                              
//! block-Jacobi-Davidson QR method. <br>                                                       
//! The generalized eigenvalues $\lambda_i$ are the diagonal entries of the                     
//! partial Schur form $A*Q = B*Q*R$ returned. <br>                                             
//!                                                                                             
//! solver options are passed in via a jadaOpts_t struct, default parameters can                
//! be set using jadaOpts_setDefaults().                                                        
//!                                                                                             
//! Input arguments:                                                                            
//!                                                                                             
//! A_op:             pointer to the operator A                                                 
//! B_op:             pointer to the hpd. operator B (if B==NULL, B=I is assumed)               
//!                                                                                             
//! for documentation of the other input parameters, see file phist_jadaOpts.h                  
//!                                                                                             
//! Output arguments:                                                                           
//!                                                                                             
//! nConv:    number of converged eigenpairs (e.g. dimension of (Q,R))                          
//! Q:        orthogonal vectors of the partial Schur form (Q,R). Should be allocated           
//!           with at least m=opts.numEigs+opts.blockSize-1 and at most opts.minBas vectors.    
//!           On output, the first nConv columns are the converged invariant subspace, but      
//!           more useful vectors (approximating yet unconverged eigenvectors) may be           
//!           contained if more vectors were allocated. A space Q of dimension opts.minBas      
//!           can be used as a starting basis for a subsequent call to subspacejada which       
//!           then does not require Arnoldi steps for starting up (set opts.v0=Q for this).     
//! R:        small upper triangular matrix of the partial Schur form (Q,R), should have the    
//!           same number of rows and columns as the number of columns in Q on input.           
//!           The upper left nConv x nConv block of R contains the Schur-form for the converged 
//!           eigenpairs, unconverged eigenvalues may be contained in the complete R if more    
//!           space is available.                                                               
//! nIter:    number of iterations performed                                                    
//! ev:       array of eigenvalues in complex arithmetic. Must be allocated by the user with    
//!           at least opts.numEigs + opts.blockSize - 1 elements.                              
//! resNorm:  norm of the residuals of the Schur form $A*q_i-Q*r_i, i=1,nEig$. Must be alloca-  
//!           ted by te user with at least opts.numEigs + opts.blockSize - 1 elements.          
//! iflag:     return code of the solver (0 on success, negative on error, positive on warning) 
//!                                                                                             
void SUBR(subspacejada)( TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) B_op,
                         phist_jadaOpts opts,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _CT_* ev,                 _MT_* resNorm,
                         int *nConv,                int *nIter,
                        int* iflag);

//! Subspace Jacobi-Davidson for interior eigenvalues

//! This algorithm is very similar to subspacejada, the main differences are    
//! that eigenvalues near a specified target inside the spectrum are sought.    
//! Instead of starting with some Arnoldi steps, we keep the shift fixed to the 
//! target until the minimum basis size is reached. The scheme computes V and W 
//! such that V'V=W'W=I, and W is an orthogonal basis of A*V (A*V=W*H_A). Eigen-
//! values are approximated using harmonic Ritz values
void SUBR(harmonicjada)( TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) B_op,
                         phist_jadaOpts opts,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _CT_ ev[],                 _MT_ resNorm[],
                         int *nConv,                int *nIter,
                        int* iflag);



//@}

//@}
