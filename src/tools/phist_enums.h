#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

typedef enum eigSort_t 
{
  NONE=0,
  LM=1,   // largest magnitude
  SM=2,   // smallest magnitude
  LR=3,   // largest real part
  SR=4,   // smallest real part
  TARGET=5, // sort according to distance from target (for interior eigenvalues)
  INVALID_EIGSORT_T=99 // returned if str2eigSort gets an invalid string
} eigSort_t;

//! how to approximately solve linear systems AX=B
typedef enum linSolv_t 
{
  DO_NOTHING=0, // do nothing: X=B
  GMRES=1, // unpreconditioned GMRES
  MINRES=2, // unpreconditioned GMRES
  CARP_CG=3, // CG on the normal equations, preconditioned by CARP (parallel SSOR)
  INVALID_LINSOLV_T=99 // returned if str2linSolv gets an invalid string
} linSolv_t;

#ifdef __cplusplus
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(eigSort_t s);
const char* linSolv2str(linSolv_t s);
eigSort_t str2eigSort(const char* str);
linSolv_t str2linSolv(const char* str);
#ifdef __cplusplus
}

//! read enum type from file stream
std::istream& operator>>(std::istream& is, eigSort_t& s);

//! read enum type from file stream
std::istream& operator>>(std::istream& is, linSolv_t& s);

#endif
#endif
