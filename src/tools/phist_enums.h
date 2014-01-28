#ifndef PHIST_ENUMS_H
#define PHIST_ENUMS_H

typedef enum eigSort_t 
{
  NONE=0,
  LM=1,   // largest magnitude
  SM=2,   // smallest magnitude
  LR=3,   // largest real part
  SR=4,   // smallest real part
  TARGET=5 // sort according to distance from target (for interior eigenvalues)
} eigSort_t;

//! how to approximately solve linear systems AX=B
typedef enum linSolv_t 
{
  DO_NOTHING=0, // do nothing: X=B
  GMRES=1, // unpreconditioned GMRES
  CARP_CG=2 // CG on the normal equations, preconditioned by CARP (parallel SSOR)
} linSolv_t;

#ifdef __cplusplus
extern "C" {
#endif
// defined in phist_tools.c
const char* eigSort2str(eigSort_t s);
const char* linSolv2str(linSolv_t s);
#ifdef __cplusplus
}
#endif
#endif
