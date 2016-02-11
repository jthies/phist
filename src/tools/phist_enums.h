#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

#include "phist_config.h"

//! defines which eigenvalues are sought resp. how they should be sorted
//! for computing eigenvalues near (far away from) a specific target, 
//! pass a shifted operator and "SM" ("LM") to the eigensolver and shift the
//! resulting eigenvalues back (eigenvectors are invariant under shifting).
typedef enum eigSort_t 
{
  phist_NO_EIGSORT=0,
  phist_LM=1,   // largest magnitude
  phist_SM=2,   // smallest magnitude
  phist_LR=3,   // largest real part
  phist_SR=4,   // smallest real part
  phist_TARGET=5, // seek eigenvalues near a specific target
  phist_INVALID_EIGSORT_T // returned if str2eigSort gets an invalid string
} eigSort_t;

//! how the eigenvalues of A are approximated in the subspace V
typedef enum eigExtr_t {
  phist_STANDARD=0, //! use Ritz values, eig(V'AV)
  phist_HARMONIC=1, //! use harmonic Ritz values, more suitable for inner eigenvalues [eig((AV)'AV)]
  phist_INVALID_EIGEXTR_T
  } eigExtr_t;

//! how to approximately solve linear systems AX=B
typedef enum linSolv_t 
{
  phist_NO_LINSOLV=0, // do nothing: X=B
  phist_GMRES=1, // unpreconditioned GMRES
  phist_MINRES=2, // unpreconditioned MINRES
  phist_CARP_CG=3, // CG on the normal equations, preconditioned by CARP (parallel SSOR)
  phist_USER_DEFINED=98,// user wants to provide custom solver
  phist_INVALID_LINSOLV_T // returned if str2linSolv gets an invalid string
} linSolv_t;

typedef enum matSym_t {
  phist_GENERAL=0,    /*! no known symmetry properties */
  phist_HERMITIAN=1,  /*! A=A^H */
  phist_COMPLEX_SYMMETRIC=2, /*! A=A^T */
  phist_PATTERN_SYMMETRIC=3, /*! G=G^T with G_ij=1 if A_ij!=0, G_ij=0 otherwise */ 
  phist_INVALID_MATSYM_T
} matSym_t;

typedef enum {
  phist_NO_PRECON=0,
  phist_IFPACK,
  phist_ML,
  phist_MUELU,
  phist_AMESOS2,
  phist_INVALID_PRECON_T
} precon_t;



#ifdef __cplusplus
#include <iostream>
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(eigSort_t s);
const char* linSolv2str(linSolv_t s);
const char* eigExtr2str(eigExtr_t s);
const char* precon2str(precon_t s);

eigSort_t str2eigSort(const char* str);
linSolv_t str2linSolv(const char* str);
eigExtr_t str2eigExtr(const char* str);
precon_t str2precon(const char* str);
#ifdef __cplusplus
}

//! read enum type from file stream
std::istream& operator>>(std::istream& is, eigSort_t& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, linSolv_t& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, eigExtr_t& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, precon_t& s);

#endif
#endif
