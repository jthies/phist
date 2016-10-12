#ifndef PHIST_JADAOPTS_H
#define PHIST_JADAOPTS_H

#include "phist_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

/*! This struct can be used to consistently pass 
    parameters to our various Jacobi-Davidson methods.
*/
typedef struct phist_jadaOpts {

// what do you want to compute?
int numEigs; //! howmany eigenpairs are sought?
phist_EeigSort which; //! LM, SM, LR, SR, or TARGET
double convTol; //! convergence tolerance for eigenvalues
phist_EmatSym symmetry; //! Symmetry properties of the matrix
phist_EeigExtr how; //! use standaard or harmonic Ritz values, etc.
               //! Generally, one should use STANDARD for extreme
               //! eigenvalues (at the border of the spectrum), and
               //! HARMONIC for inner ones. Other methods may be
               //! added later.

/**********************************
 * JaDa configuration             *
 **********************************/

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

/**********************************
 * inner solver configuration     *
 **********************************/

phist_ElinSolv innerSolvType; //! GMRES, MINRES, CARP_CG, USER_DEFINED currently supported.
                              //! If set to USER_DEFINED, you have to provide the customSolver*
                              //! interface below.

int innerSolvBlockSize;
int innerSolvMaxBas;
int innerSolvMaxIters;
int innerSolvRobust; /*! extra effort to get good jada updates
                      * (in practice this may mean a more accurate orthogonalization etc.)
                      */
  int innerSolvStopAfterFirstConverged;
  
  int innerSolvMaxProjectionSpace; /* allow at most <k> vectors to be projected out in the
                                    * inner solver, even if more eigenvectors have been locked already.
                                    A value <0 indicates that all locked vectors are projected out.
                                    */

  //! pointer to a inearOp whose apply_shifted function serves as a preconditioner for the inner solver. May be NULL (no preconditioning) or created by phist_Xprecon_create.
  //! This pointer can be used to pass an already computed preconditioner to the Jacobi-Davidson solver. The preconditioner will not be updated explicitly during the run, but
  //! differentsshift will be supplied to the apply_shifted function.
  void const* preconOp;

  //! This fields allows specifying a preconditioner (along with preconOpts to pass in options) in an option file.
  //! Currently it is the responsibility of the driver routine to create the preconditioner in "preconOp", but in 
  //! the future we may allow the jada solver(s) to update the preconditioner with new subspace information and   
  //! shifts during the eigenvalue computation.
  phist_Eprecon preconType;

  //! option string passed to precon_create alongside preconType (if it is not NO_PRECON or INVALID_PRECON)
  char preconOpts[1024];

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

//! This function allows customizing the residual calculation jdqr, it is a tailored
//! interface for the HYMLS (http://bitbucket.org/hymls/hymls) project and not intended
//! for other uses right now. Setting this does *not* affect subspacejada.
void (*custom_computeResidual)(void* customSolverData, void const* B_op, void* r_ptr,
        void const* Au_ptr, void const* u_ptr, void* rtil_ptr,
        void const* Qv, void* tmp, void* Theta,
        void* atil, void* *atilv, double *resid,
        int nv, int nconv, int* iflag);


} phist_jadaOpts;

//! get jada options from a simple ASCII file and place them in the struct passed to the solvers.

//! overwrite the struct members from a file containing lines like this:
//!
//! numEigs 8
//! which LM
//! how STANDARD
//! convTol 1.0e-8
//!
//! The function is not at all fancy, it won't warn about invalid entries, doesn't
//! care about invalid lines like BLABLA_numEigs -99 and may throw exceptions, especially
//! if the file can't be opened. Every MPI process opens the file separately.
void phist_jadaOpts_fromFile(phist_jadaOpts *opts, const char* filename, int* iflag);

//! set default values in a jadaOpts struct
void phist_jadaOpts_setDefaults(phist_jadaOpts *opts);

//! print jadaOpts to a file or stream. The result can be used as input for subsequent runs
void phist_jadaOpts_toFile(phist_jadaOpts const *opts, FILE* stream);

#ifdef __cplusplus
} //extern "C"
#endif

#endif
