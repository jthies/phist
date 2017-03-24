/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! internal implementation of the orthog function, the regular SUBR(orthog) should be used instead in actual algorithms.

//! Does not randomize the null space (if any). Note that the arguments are renamed and reordered, this should 
//! eventually be adjusted to avoid confusion. Here, W is the orthogonal space, V is orthogonalized against W, and the
//! relation that holds after the subroutine (with Q=V on output) is Q*R1 = V-W*R2, Q'Q=I.
//!
//! The return value of this function in *iflag is the rank of the null space of [W V] on input. The last *iflag
//! columns of V are 0.
void SUBR(orthogrrfused)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, 
        TYPE(const_sdMat_ptr) VtV, int* iflag);
