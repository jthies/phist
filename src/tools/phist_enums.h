#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

//! defines which eigenvalues are sought resp. how they should be sorted
//! for computing eigenvalues near (far away from) a specific target, 
//! pass a shifted operator and "SM" ("LM") to the eigensolver and shift the
//! resulting eigenvalues back (eigenvectors are invariant under shifting).
typedef enum eigSort_t 
{
  NONE=0,
  LM=1,   // largest magnitude
  SM=2,   // smallest magnitude
  LR=3,   // largest real part
  SR=4,   // smallest real part
  TARGET=5, // seek eigenvalues near a specific target
  INVALID_EIGSORT_T=99 // returned if str2eigSort gets an invalid string
} eigSort_t;

//! how the eigenvalues of A are approximated in the subspace V
typedef enum eigExtr_t {
  STANDARD=0, //! use Ritz values, eig(V'AV)
  HARMONIC=1, //! use harmonic Ritz values, more suitable for inner eigenvalues [eig((AV)'AV)]
  INVALID_EIGEXTR_T=99
  } eigExtr_t;

//! how to approximately solve linear systems AX=B
typedef enum linSolv_t 
{
  DO_NOTHING=0, // do nothing: X=B
  GMRES=1, // unpreconditioned GMRES
  MINRES=2, // unpreconditioned MINRES
  CARP_CG=3, // CG on the normal equations, preconditioned by CARP (parallel SSOR)
  USER_DEFINED=98,// user wants to provide custom solver
  INVALID_LINSOLV_T=99 // returned if str2linSolv gets an invalid string
} linSolv_t;

typedef enum matSym_t {
  GENERAL=0,    /*! no known symmetry properties */
  HERMITIAN=1,  /*! A=A^H */
  COMPLEX_SYMMETRIC=2, /*! A=A^T */
  PATTERN_SYMMETRIC=3, /*! G=G^T with G_ij=1 if A_ij!=0, G_ij=0 otherwise */ 
  INVALID_MATSYM_T=99
} matSym_t;

#ifdef __cplusplus
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(eigSort_t s);
const char* linSolv2str(linSolv_t s);
const char* eigExtr2str(eigExtr_t s);
eigSort_t str2eigSort(const char* str);
linSolv_t str2linSolv(const char* str);
eigExtr_t str2eigExtr(const char* str);
#ifdef __cplusplus
}

//! read enum type from file stream
std::istream& operator>>(std::istream& is, eigSort_t& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, linSolv_t& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, eigExtr_t& s);

#endif
#endif
