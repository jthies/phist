/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_jadaCorrectionSolver_decl.h 
//! \brief calculate approximate solutions to a set of Jacobi-Davidson correction equations

//! \addtogroup jada
//!@{

//! \brief The jadaCorrectionSolver uses the blockedGMRES to calculate approximate solutions to a set of Jacobi-Davidson correction equations.
//!
//! It provides a simple interface and takes care of restarting/pipelining issues
//!
typedef struct TYPE(jadaCorrectionSolver) {
  //! \name Internal data structures
  //!@{
  int                   hermitian; //!< Some inner solver schemes may be able to make use of this info
  int                   innerSolvBlockSize_; //!< Number of blockedGMRES states iterated at once
  int                   innerSolvMaxBas_; //!< number of blocks allowed in linear solver (m for GMRES(m), s for IDR(s))
  TYPE(blockedGMRESstate_ptr) *blockedGMRESstates_; //!< blockedGMRES states
  TYPE(carp_cgState_ptr) *carp_cgStates_; //!< Can use CARP-CG alternatively
  phist_ElinSolv     method_;    //!< Supported values are GMRES, MINRES, CARP_CG and CUSTOM.
  
  int preconSkewProject;

  TYPE(linearOp_ptr) leftPrecon;

  TYPE(linearOp_ptr) rightPrecon;

  //! Pointer to solver object if innerSolvType==USER_DEFINED
  void* customSolver_;

  //! \brief This function is used instead of phist_jadaCorrectionSolver_run if innerSolvType is USER_DEFINED.
  //!
  //! For subspacejada with block size 1 it is enough to implement the simpler interface below,
  //! we first check in those cases if that interface is set before checking for this one.
  void (*customSolver_run)(         void*  customSolverData,
                                    void const*    A_op,       void const*    B_op,
                                    void const*    Qtil,       void const*    BQtil,
                                    const double sigma_r[],    const double sigma_i[],  
                                    TYPE(const_mvec_ptr)  res, const int resIndex[],
                                    const double        tol[], int       maxIter,
                                    TYPE(mvec_ptr)        t,
                                    int robust,                int abortAfterFirstConvergedInBlock,
                                    int *iflag);

  //! Simplified interface if only single-vector jdqr or subspacejada is used.
  void (*customSolver_run1)(        void*  customSolverData,
                                    void const*    A_op,     void const*    B_op,
                                    void const*    Qtil,     void const*    BQtil,
                                    double   sigma_r,        double sigma_i, 
                                     TYPE(const_mvec_ptr)  res,
                                    const double           tol,    int                   maxIter,
                                    TYPE(mvec_ptr)        t,
                                    int robust,
                                    int *                 iflag);
  //!@}
} TYPE(jadaCorrectionSolver);

//! Pointer to jadaCorrectionSolver object
typedef TYPE(jadaCorrectionSolver)* TYPE(jadaCorrectionSolver_ptr);
//! Pointer to const jadaCorrectionSolver object
typedef TYPE(jadaCorrectionSolver) const * TYPE(const_jadaCorrectionSolver_ptr);


//! Create a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_create)(TYPE(jadaCorrectionSolver_ptr) *me, phist_jadaOpts opts,
        phist_const_map_ptr map, int *iflag);
        
//! Delete a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_delete)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver, int *iflag);


//! \brief Calculate approximate solutions to given set of Jacobi-Davidson correction equations
//!
//! \param jdCorrSolver    the jadaCorrectionSolver object
//! \param AB_op            matrix A passed to jadaOp_create
//! \param B_op            matrix B passed to jadaOp_create
//! \param Qtil            projection vectors V passed to jadaOp_create
//! \param BQtil           projection vectors BV passed to jadaOp_create
//! \param sigma           (pos.!) shifts, -sigma[i], i in {1, ..., nvec} is passed to the jadaOp
//! \param res             JD residuals, e.g. rhs of the correction equations
//! \param resIndex        if not NULL, specifies permutation of the residual array to avoid unnecessary copying in the jada-algorithm
//! \param tol             desired accuracy (gmres residual tolerance) of the individual systems
//! \param maxIter         maximal number of iterations after which individial systems should be aborted
//! \param t               returns approximate solution vectors
//! \param iflag            a value > 0 indicates the number of systems that have not converged to the desired tolerance
void SUBR(jadaCorrectionSolver_run)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver,
                                    TYPE(const_linearOp_ptr)    AB_op, TYPE(const_linearOp_ptr) B_op, 
                                    TYPE(const_mvec_ptr)  Qtil,        TYPE(const_mvec_ptr)     BQtil,
                                    const _ST_            sigma[],     TYPE(const_mvec_ptr)     res,      
                                    const int             resIndex[],  const _MT_               tol[],    
                                    int                   maxIter,     TYPE(mvec_ptr)           t,
                                    int useIMGS, int abortAfterFirstConvergedInBlock, int updatePrecon,
                                    int *                 iflag);

//!@}
