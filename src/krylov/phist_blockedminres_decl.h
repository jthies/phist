/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \addtogroup linear_solvers
//@{

//! \defgroup minres MinRes solver for symmetric indefinite matrices
//@{

//! perform minres iterations (for symmetric indefinite matrices)

//! similar to blockedGMRESstates_iterate, but uses MINRES (e.g. the Lanczos procedure instead of Arnoldi)
//! suited only for symmetric/Hermitian matrices
void SUBR( blockedMINRESstates_iterate ) (TYPE(const_linearOp_ptr) Op, 
                                          TYPE(const_linearOp_ptr) rightPrecon,
                                          TYPE(blockedGMRESstate_ptr) S_array[], int numSys, int *nIter, int* iflag);

//@}

//@}
