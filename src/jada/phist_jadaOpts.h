/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_JADAOPTS_H
#define PHIST_JADAOPTS_H

#ifndef DOXYGEN
#include "phist_enums.h"
#endif

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
int lookAhead; //! consider at most 'lookAhead' unconverged eigenvalues at a time (set it to -1 for the default of 2*blockSize)
int minBas; //! number of vectors retained upon restart
int maxBas; //! maximum number of vectors allowed in the basis

// how should JaDa start up?
void* v0; //! can be used to pass in a start-up vector(-space) (can have any number of 
          //! columns). v0 is assumed to be orthonormal. Our current implementation of 
          //! subspacejada will use only the first column of v0 and perform Arnoldi iterations
          //! unless v0 has exactly minBas columns. In this case v0 is used instead of
          //! running Arnoldi iterations. A suitable input basis can be obtained from a previous
          //! run of subspacejada by pre-allocating Q with minBas vectors.
int arno; //! 0: no Arnoldi steps. 1: minBas Arnoldi steps to start up.
double initialShift_r; //! can be used to start with an initial shift
                       //! (ignored if arno!=0)
double initialShift_i; //! imaginary part of initial shift

int initialShiftIters; // perform given number of iterations with a fixed shift

/**********************************
 * inner solver configuration     *
 **********************************/

phist_ElinSolv innerSolvType; //! GMRES, MINRES, CARP_CG, USER_DEFINED and NO_LINSOLV currently supported.
                              //! If set to USER_DEFINED, you have to provide the customSolver*
                              //! interface below.
                              //! The value "NO_LINSOLV", combined with a preconditioner, leads to an Olsen
                              //! method (only the preconditioner is applied, with appropriate projections 
                              //! if preconSkewProject!=0 is set.

int innerSolvBlockSize;       //! if set to -1, the outer block size (blockSize) is used.
int innerSolvMaxBas;          //! if set to -1, innerSolvMaxBas=innerSolvMaxIters is used (no restarting of e.g. GMRES)
int innerSolvMaxIters;
double innerSolvBaseTol;    //! in the k'th iteration on eigenvalue j, the tolerance used will be innerTolBase^k (default: 0.1)
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
  //!
  //! note: this field has no 'const' qualifier, however, if preconUpdate==0, the preconditioner pointed to
  //! is in fact not modified.
  void* preconOp;

  //! This field allows specifying a preconditioner (along with preconOpts to pass in options) in an option file.
  //! It is the responsibility of the driver routine to create the preconditioner in "preconOp"
  phist_Eprecon preconType;

  //! option string passed to precon_create alongside preconType (if it is not NO_PRECON or INVALID_PRECON)
  char preconOpts[1024];
  
  //! if 0, just apply the preconditioner K "as is"
  //! if 1, apply skew-projection with the currend approximation(s) q, giving the operator
  //!       (I-(K\q) ((Bq)'(K\q))^{-1}(Bq)')K^{-1}.
  //! if 2, a symmetric variant with pre- and postprojection). This is currently not imple-
  //!       mented, though, because we currently only have MINRES to exploit symmetry in the
  //!       inner iterations. Unless all locked eigenvectors are projected out, the resulting
  //!       operator will be indefinite even if K is hpd, so it can't be used as a precon-
  //!       ditioner for MINRES anyway.
  int preconSkewProject;
  
  //! if 0, the preconditioner is kept the same throughout the
  //! Jacobi-Davidson process.
  //! If 1 it is updated once with as shift the initial approximate eigenvalue obtained
  //! by the starting space v0 or Arnoldi iterations, respectively.
  //! If >1, it is updated with the current shift before each correction solve.
  //!
  //! Variants like updating based on convergence rate may be added later.
  int preconUpdate;

  //! pointer to solver object if innerSolvType==USER_DEFINED
  void* customSolver;

  //! this function is used instead of phist_jadaCorrectionSolver_run if innerSolvType is USER_DEFINED.
  //! For subspacejada with block size 1 it is enough to implement the simpler interface below,
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

  //! simplified interface if only single-vector subspacejada is used.
  void (*customSolver_run1)(        void*  customSolverData,
                                    void const* A_op,     void const*    B_op,
                                    void const* Qtil,     void const*    BQtil,
                                    const double sigma,   double sigma_i,
                                    void const* res,      double tol,
                                    int maxIter,          void* t,
                                    int robust,
                                    int * iflag);

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
