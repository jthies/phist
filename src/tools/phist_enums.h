/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

#include "phist_config.h"

//! defines which eigenvalues are sought resp. how they should be sorted
//! for computing eigenvalues near (far away from) a specific target, 
//! pass a shifted operator and "SM" ("LM") to the eigensolver and shift the
//! resulting eigenvalues back (eigenvectors are invariant under shifting).
typedef enum phist_EeigSort 
{
  phist_NO_EIGSORT=0,
  phist_LM=1,   // largest magnitude
  phist_SM=2,   // smallest magnitude
  phist_LR=3,   // largest real part
  phist_SR=4,   // smallest real part
  phist_TARGET=5, // seek eigenvalues near a specific target
  phist_INVALID_EIGSORT // returned if str2eigSort gets an invalid string
} phist_EeigSort;

//! how the eigenvalues of A are approximated in the subspace V
typedef enum phist_EeigExtr 
{
  phist_STANDARD=0, //! use Ritz values, eig(V'AV)
  phist_HARMONIC=1, //! use harmonic Ritz values, more suitable for inner eigenvalues [eig((AV)'AV)]
  phist_INVALID_EIGEXTR
} phist_EeigExtr;

//! how to approximately solve linear systems AX=B
typedef enum phist_ElinSolv 
{
  phist_NO_LINSOLV=0, // do nothing: X=B
  phist_GMRES=1, // GMRES for general matrices (with left and/or right preconditioning)
  phist_MINRES=2, // MINRES for symmetric and indefinite matrices
  phist_QMR=3, // QMR for symmetric and indefinite matrices
  phist_CARP_CG=4, // CG on the normal equations, preconditioned by CARP (parallel SSOR)
  phist_USER_LINSOLV=98,// user wants to provide custom solver
  phist_INVALID_LINSOLV // returned if str2linSolv gets an invalid string
} phist_ElinSolv;

typedef enum phist_EmatSym 
{
  phist_GENERAL=0,    /*! no known symmetry properties */
  phist_HERMITIAN=1,  /*! A=A^H */
  phist_COMPLEX_SYMMETRIC=2, /*! A=A^T */
  phist_PATTERN_SYMMETRIC=3, /*! G=G^T with G_ij=1 if A_ij!=0, G_ij=0 otherwise */ 
  phist_INVALID_MATSYM
} phist_EmatSym;

typedef enum 
{
  phist_NO_PRECON=0,
  phist_IFPACK,
  phist_ML,
  phist_MUELU,
  phist_AMESOS2,
  phist_USER_PRECON=98,
  phist_INVALID_PRECON
} phist_Eprecon;



#ifdef __cplusplus
#include <iostream>
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(phist_EeigSort s);
const char* linSolv2str(phist_ElinSolv s);
const char* eigExtr2str(phist_EeigExtr s);
const char* precon2str(phist_Eprecon s);
const char* matSym2str(phist_EmatSym s);

phist_EeigSort str2eigSort(const char* str);
phist_ElinSolv str2linSolv(const char* str);
phist_EeigExtr str2eigExtr(const char* str);
phist_Eprecon str2precon(const char* str);
phist_EmatSym str2matSym(const char* str);
#ifdef __cplusplus
}

//! read enum type from file stream
std::istream& operator>>(std::istream& is, phist_EeigSort& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, phist_ElinSolv& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, phist_EeigExtr& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, phist_Eprecon& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, phist_EmatSym& s);

#endif
#endif
