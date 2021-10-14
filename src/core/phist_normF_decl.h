/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_normF_decl.h
//! \brief Implementation of Fobenius norm

//! \name Frobenius norm
//! \addtogroup core
//!@{

//!\brief Frobenius norm for sdMats
//!
//! For GPU processes, will run on the host side only, so the caller may have to
//! sync sdMat data using sdMat_from_device.
void SUBR(sdMat_normF)(TYPE(const_sdMat_ptr) M, _MT_ *f, int* iflag);

//!\brief Frobenius norm for mvecs
void SUBR(mvec_normF)(TYPE(const_mvec_ptr) V, _MT_ *f, int* iflag);

//!@}

