/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_kpm_decl.h
//! \brief Implementation of KPM algorithms

//! \addtogroup poly
//!@{

//! \brief KPM algorithm (Weisse, Wellein, Alvermann & Fehske, Rev. Mod. Phys. 78 (1), 275-306).
//!
//! \param scale,shift have to map the spectrum of A into [-1,1].
//! \param M is the number of iterations.
//! \param [out] mu is assumed to have size 2 * M * |cols(x)|.
//!
//! Currently there are no additional options for iflag.
void SUBR(kpm)(TYPE(mvec_ptr) x, TYPE(sparseMat_ptr) A, _MT_ scale, _MT_ shift, int M, _MT_* mu, int* iflag);

//! \brief KPM-DOS algorithm (Weisse, Wellein, Alvermann & Fehske, Rev. Mod. Phys. 78 (1), 275-306).
//!
//! \param scale,shift have to map the spectrum of A into [-1,1].
//! \param M is the number of iterations.
//! \param R is the number of vectors.
//! \param [out] mu is assumed to have size 2 * M.
//!
//! The flag PHIST_KPM_SINGLEVEC iterates a single vector multiple times instead of computing multiple vectors simultaneously.
void SUBR(dos)(TYPE(sparseMat_ptr) A, _MT_ scale, _MT_ shift, int M, int R, _MT_* mu, int* iflag);

//!@}
