/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_CORE_FLAGS_H
#define PHIST_CORE_FLAGS_H

/*! \file phist_core_flags.h
   
   some flags that influence the behavior of core functions, see also
   phist_kernel_flags.h for an introduction to flags in PHIST.
   
   The values 2^[16...23] are reserved for core functionality for now, again, 
   if there is no function that accepts flags A and B, A==B is allowed.
   
 */

 //! \name some flags that influence the behavior of core functions
 //! \addtogroup core
 //!@{
 
//! \def PHIST_KPM_SINGLEVEC
//!
//!    Accepted by KPM core routine 
#define PHIST_KPM_SINGLEVEC 65536

//! \def PHIST_ORTHOG_RANDOMIZE_NULLSPACE
//!  
//!  For orthog routines: fill up output vector with random numbers and orthogonalize them along 
//!   with the original entries if the vector W is found to be rank-deficient.                    
#define PHIST_ORTHOG_RANDOMIZE_NULLSPACE 131072

//! \def PHIST_ORTHOG_TRIANGULAR_R1
//!   
//!    For orthog routine: make sure the resulting factor R1 is upper triangular on output,
//!    orthog computes [Q,R1,R2] s.t. Q*R1 = W - V*R2, with V an orthogonal basis on input.
#define PHIST_ORTHOG_TRIANGULAR_R1       262144

 //!@}

#endif
