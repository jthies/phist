/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_subspacejada_decl.h
//! \brief Jacobi-Davidson QR method with blocking, locking and restart

//! \brief Subspace Jacobi-Davidson for extreme eigenvalues \ingroup bjdqr

//! Tries to compute a partial Schur form \f$(Q,R)\f$ of dimension opts.nEigs                       
//! of the stencil \f$A*x-\lambda*B*x\f$ with a general linear operator \f$A\f$ and a                   
//! hermitian positive definite (hpd.) linear operator \f$B\f$ using a                              
//! block-Jacobi-Davidson QR method. <br>                                                       
//! The generalized eigenvalues \f$\lambda_i\f$ are the diagonal entries of the                     
//! partial Schur form \f$A*Q = B*Q*R\f$ returned. <br>                                             
//!                                                                                             
//! Solver options are passed in via a jadaOpts_t struct, default parameters can                
//! be set using jadaOpts_setDefaults().                                                        
//!                                                                                                                                                                                          
//! \param [in] A_op  pointer to the operator A                                                 
//! \param [in] B_op  pointer to the hpd. operator B (if B==NULL, B=I is assumed)               
//!                                                                                             
//! For documentation of the other input parameters, see file phist_jadaOpts.h                                                                                            
//!                                                                                             
//! \param [out] nConv    number of converged eigenpairs (e.g. dimension of (Q,R))                          
//! \param [out] Q        orthogonal vectors of the partial Schur form (Q,R). Should be allocated           
//!                       with at least m=opts.numEigs+opts.blockSize-1 and at most opts.minBas vectors. <br>   
//!                       On output, the first nConv columns are the converged invariant subspace, but      
//!                       more useful vectors (approximating yet unconverged eigenvectors) may be           
//!                       contained if more vectors were allocated. A space Q of dimension opts.minBas      
//!                       can be used as a starting basis for a subsequent call to subspacejada which       
//!                       then does not require Arnoldi steps for starting up (set opts.v0=Q for this). <br>    
//!                       It is allowed to use the same Q as v0 (input) and output argument.                
//! \param [out] R        small upper triangular matrix of the partial Schur form (Q,R), should have the    
//!                       same number of rows and columns as the number of columns in Q on input. <br>
//!                       The upper left nConv x nConv block of R contains the Schur-form for the converged 
//!                       eigenpairs, unconverged eigenvalues may be contained in the complete R if more    
//!                       space is available.                                                               
//! \param [out] nIter    number of iterations performed                                                    
//! \param [out] ev       array of eigenvalues in complex arithmetic. Must be allocated by the user with    
//!                       at least opts.numEigs + opts.blockSize - 1 elements.                              
//! \param [out] resNorm  norm of the residuals of the Schur form \f$A*q_i-Q*r_i, i=1,nEig\f$. Must be alloca-  
//!                       ted by te user with at least opts.numEigs + opts.blockSize - 1 elements.          
//! \param [out] iflag    return code of the solver (0 on success, negative on error, positive on warning) 
//!   
//! For implementation of the function, see file src/jada/phist_subspacejada_def.hpp
void SUBR(subspacejada)( TYPE(const_linearOp_ptr) A_op,  TYPE(const_linearOp_ptr) B_op,
                         phist_jadaOpts opts,
                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
                         _CT_* ev,                 _MT_* resNorm,
                         int *nConv,                int *nIter,
                        int* iflag);



