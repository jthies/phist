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
double convTol; //! convergence tolerance for eigenvalues
matSym_t symmetry; //! Symmetry properties of the matrix

////////////////////////////////////
// JaDa configuration             //
////////////////////////////////////

int maxIters; //! maximum iterations allowed
int blockSize; //! only for block methods (subspacejada)
int minBas; //! number of vectors retained upon restart
int maxBas; //! maximum number of vectors allowed in the basis

// how should JaDa start up?
void* v0; //! can be used to pass in a start-up vector(-space) (can have any number of 
          //! columns). v0 is assumed to be orthonormal.
int arno; //! 0: no Arnoldi steps. 1: minBas Arnoldi steps to start up.
double initialShift_r; //! can be used to start with an initial shift
                       //! (ignored if arno!=0)
double initialShift_i; //! imaginary part of initial shift

int initialShiftIters; // perform given number of iterations with a fixed shift

/*//////////////////////////////////
// inner solver configuration     //
//////////////////////////////////*/

linSolv_t innerSolvType; /*! GMRES, MINRES, CARP_CG, USER_DEFINED currently supported.
                          * If set to USER_DEFINED, you have to provide the customSolver*
                          * interface below.
                          */
int innerSolvBlockSize;
int innerSolvMaxBas;
int innerSolvMaxIters;
int innerSolvRobust; /*! extra effort to get good jada updates
                      * (in practice this may mean a more accurate orthogonalization etc.)
                      */
int innerSolvStopAfterFirstConverged;


  //! pointer to solver object if innerSolvType==USER_DEFINED
  void* customSolver;

  //! this function is used instead of phist_jadaCorrectionSolver_run if innerSolvType is USER_DEFINED.
  //! For jdqr or subspacejada with block size 1 it is enough to implement the simpler interface below,
  //! we first check in those cases if that interface is set before checking for this one.
  //! note that the scalar arguments are always passed in as doubles so that a single function pointer can
  //! be used in this untyped struct.
  void (*customSolver_run)(         void*       customSolverData,
                                    void const* A_op,      void const*    B_op,
                                    void const* Qtil,      void const*    BQtil,
                                    const double sigma_r[], const double sigma_i[],  
                                    void const* res,        const int resIndex[],
                                    const double tol[],       int maxIter,
                                    void* t,                int robust,
                                    int abortAfterFirstConvergedInBlock,
                                    int * iflag);

  //! simplified interface if only single-vector jdqr or subspacejada is used.
  void (*customSolver_run1)(        void*  customSolverData,
                                    void const* A_op,     void const*    B_op,
                                    void const* Qtil,     void const*    BQtil,
                                    const double sigma,   double sigma_i,
                                    void const* res,      double tol,
                                    int maxIter,          void* t,
                                    int robust,
                                    int * iflag);

} phist_jadaOpts_t;

void phist_jadaOpts_setDefaults(phist_jadaOpts_t *opts);

#ifdef __cplusplus
} //extern "C"
#endif

#endif
