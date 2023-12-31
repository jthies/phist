/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_enums.h
//! \brief phist enums 
#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

#include "phist_config.h"

//! \defgroup enums phist enums
//! \ingroup tools
//!@{

//! \brief defines which eigenvalues are sought resp. how they should be sorted
//!
//! for computing eigenvalues near (far away from) a specific target, 
//! pass a shifted operator and "SM" ("LM") to the eigensolver and shift the
//! resulting eigenvalues back (eigenvectors are invariant under shifting).
typedef enum phist_EeigSort 
{
  phist_NO_EIGSORT=0,
  phist_LM=1,   //!< largest magnitude
  phist_SM=2,   //!< smallest magnitude
  phist_LR=3,   //!< largest real part
  phist_SR=4,   //!< smallest real part
  phist_TARGET=5, //!< seek eigenvalues near a specific target
  phist_INVALID_EIGSORT //!< returned if str2eigSort gets an invalid string
} phist_EeigSort;

//! how the eigenvalues of A are approximated in the subspace V
typedef enum phist_EeigExtr 
{
  phist_STANDARD=0, //!< use Ritz values, eig(V'AV)
  phist_HARMONIC=1, //!< use harmonic Ritz values, more suitable for inner eigenvalues [eig((AV)'AV)]
  phist_INVALID_EIGEXTR
} phist_EeigExtr;

//! how to approximately solve linear systems AX=B
typedef enum phist_ElinSolv 
{
  phist_NO_LINSOLV=0, //!< do nothing: X=B
  phist_GMRES=1, //!< GMRES for general matrices (with left and/or right preconditioning)
  phist_MINRES=2, //!< MINRES for symmetric and indefinite matrices
  phist_QMR=3, //!< QMR for general matrices
  phist_BICGSTAB=4, //!< BiCGStab for general matrices
  phist_CARP_CG=5, //!< CG on the normal equations, preconditioned by CARP (parallel SSOR)
  phist_PCG=6, //!< preconditioned Conjugate Gradients for Hermitian systems
  phist_IDRS=7, //!< preconditioned IDR(s) for general systems
  phist_USER_LINSOLV=98,//!< user wants to provide custom solver
  phist_INVALID_LINSOLV //!< returned if str2linSolv gets an invalid string
} phist_ElinSolv;

//! symmetry of matrix A
typedef enum phist_EmatSym 
{
  phist_GENERAL=0,    //!< no known symmetry properties
  phist_HERMITIAN=1,  //!< A=A^H
  phist_COMPLEX_SYMMETRIC=2, //!< A=A^T
  phist_PATTERN_SYMMETRIC=3, //! G=G^T with G_ij=1 if A_ij!=0, G_ij=0 otherwise
  phist_INVALID_MATSYM
} phist_EmatSym;

//! preconditioner
typedef enum phist_Eprecon
{
  phist_NO_PRECON=0,
  phist_IFPACK,
  phist_ML,
  phist_MUELU,
  phist_AMESOS2,
  phist_USER_PRECON=98,
  phist_INVALID_PRECON
} phist_Eprecon;

//! \brief define operators involving pre, post and skew projection
//! for the Jacobi-Davidson method.
//!
//! The general form of the jadaOp (see phist_jadaOp_decl.h) is
//!
//! P_skew K^{-1} P_post (A-sigma_j B) P_pre,
//!
//! with W=B*V (B hpd, B=I for standard EVP),
//! K some preconditioner for (A-sigma_j *B),
//! V_K=K\ V, and (') denoting transposition, 
//! the operators are defined as
//!
//! P_skew = I - V_K*(W'V_K)^{-1}W' <br>
//! P_pre  = I - V (BV)' <br>
//! P_post = P_pre' = I - BV V'
//!
//! Using this enum one can 'turn on/off' the projections,
//! for insance phist_PROJ_SKEW|phist_PROJ_PRE yields the operator
//! P_skew K^{-1} (A-sigma_j B) P_pre.
typedef enum phist_Eprojection
{
  phist_PROJ_NONE=0,
  phist_PROJ_PRE=1,
  phist_PROJ_POST=2,
  phist_PROJ_PRE_POST=3,
  phist_PROJ_SKEW=4,
  phist_PROJ_ALL=7,
  phist_INVALID_PROJ
} phist_Eprojection;

#ifdef __cplusplus
# ifndef DOXYGEN
# include <iostream>
# endif
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(phist_EeigSort s);
const char* linSolv2str(phist_ElinSolv s);
const char* eigExtr2str(phist_EeigExtr s);
const char* precon2str(phist_Eprecon s);
const char* matSym2str(phist_EmatSym s);
const char* projection2str(phist_Eprojection s);

phist_EeigSort str2eigSort(const char* str);
phist_ElinSolv str2linSolv(const char* str);
phist_EeigExtr str2eigExtr(const char* str);
phist_Eprecon str2precon(const char* str);
phist_EmatSym str2matSym(const char* str);
phist_Eprojection str2projection(const char* str);

#ifdef __cplusplus
}

//! \name read enum types from input stream \ingroup enums
//! @{

//!
std::istream& operator>>(std::istream& is, phist_EeigSort& s);
//!
std::istream& operator>>(std::istream& is, phist_ElinSolv& s);
//!
std::istream& operator>>(std::istream& is, phist_EeigExtr& s);
//!
std::istream& operator>>(std::istream& is, phist_Eprecon& s);
//!
std::istream& operator>>(std::istream& is, phist_EmatSym& s);
//!
std::istream& operator>>(std::istream& is, phist_Eprojection& s);

//!@}

//! \name write enum types to stream \ingroup enums
//! @{

//!
std::ostream& operator<<(std::ostream& is, const phist_EeigSort& s);
//!
std::ostream& operator<<(std::ostream& is, const phist_ElinSolv& s);
//!
std::ostream& operator<<(std::ostream& is, const phist_EeigExtr& s);
//!
std::ostream& operator<<(std::ostream& is, const phist_Eprecon& s);
//!
std::ostream& operator<<(std::ostream& is, const phist_EmatSym& s);
//!
std::ostream& operator<<(std::ostream& is, const phist_Eprojection& s);

//!@}

#endif
//!@}
#endif
