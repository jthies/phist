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
                         _CT_* ev,                 _MT_* resNorm,
                         int *nConv,                int *nIter,
                        int* iflag);



//@}

//@}
