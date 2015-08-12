#ifndef PHIST_JADAOPTS_H
#define PHIST_JADAOPTS_H

#include "phist_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! This struct can be used to consistently pass 
    parameters to our various Jacobi-Davidson methods.
*/
typedef struct phist_jadaOpts_t {

// what do you want to compute?
int numEigs; //! howmany eigenpairs are sought?
eigSort_t which; //! LM, SM, LR, SR, or TARGET
linSolv_t innerSolvType; //! GMRES, MINRES, CARP_CG, USER_DEFINED currently supported.
                         //! If set to USER_DEFINED, you have to provide the customSolver*
                         //! interface below.
double convTol; //! convergence tolerance for eigenvalues

// JaDa configuration
int maxIters; //! maximum iterations allowed
int blockSize; //! only for block methods (subspacejada and blockjada)
int minBas; //! number of vectors retained upon restart
int maxBas; //! maximum number of vectors allowed in the basis

// how should JaDa start up?
void* v0; //! can be used to pass in a start-up vector(-space) (can have any number of 
          //! columns). v0 is assumed to be orthonormal.
int arno; //! 0: no Arnoldi steps. 1: minBas Arnoldi steps to start up.
double initialShift; //! can be used to start with an initial shift
                     //! (ignored if arno!=0)

  //! custom solver data structure if innerSolvType is USER_DEFINED
void* customSolver_;
#if 0  
  void (*customSolver_create)( void** customSolverData,int blockDim, const_map_ptr_t map, int maxBase, int* iflag);
  void (*customSolver_delete)( void*  customSolverData, int* iflag);
  void (*customSolver_run)(            void*  customSolverData,
                                    void const*    A_op,     void const*    B_op,
                                    void const*    Qtil,     void const*    BQtil,
                                    const _ST_            sigma[],  TYPE(const_mvec_ptr)  res,      const int resIndex[],
                                    const _MT_            tol[],    int                   maxIter,
                                    TYPE(mvec_ptr)        t,
                                    bool useIMGS,                   bool abortAfterFirstConvergedInBlock,
                                    int *                 iflag);
#endif  

} phist_jadaOpts_t;

void phist_jadaOpts_setDefaults(phist_jadaOpts_t *opts);

#ifdef __cplusplus
} //extern "C"
#endif

#endif
