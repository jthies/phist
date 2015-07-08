#include "phist_enums.h"

//!
//! The jadaCorrectionSolver uses the blockedGMRES to calculate approximate solutions to a set of Jacobi-Davidson correction equations.
//! It provides a simple interface and takes care of restarting/pipelining issues
//!
typedef struct TYPE(jadaCorrectionSolver)
{
  //! \name internal data structures
  //@{
  int                   gmresBlockDim_;     //! number of blockedGMRES states iterated at once
  TYPE(blockedGMRESstate_ptr) *blockedGMRESstates_;     //! blockedGMRES states
  TYPE(carp_cgState_ptr) *carp_cgStates_; //! can use CARP-CG alternatively
  linSolv_t     method_;    //! supported values are GMRES, MINRES and CARP_CG
  //@}
} TYPE(jadaCorrectionSolver);

typedef TYPE(jadaCorrectionSolver)* TYPE(jadaCorrectionSolver_ptr);

typedef TYPE(jadaCorrectionSolver) const * TYPE(const_jadaCorrectionSolver_ptr);


//! create a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_create)(TYPE(jadaCorrectionSolver_ptr) *jdCorrSolver, int blockedGMRESBlockDim, const_map_ptr_t map, 
        linSolv_t method, int maxBase, int *iflag);

//! delete a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_delete)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver, int *iflag);


//! calculate approximate solutions to given set of jacobi-davidson correction equations
//!
//! arguments:
//! jdCorrSolver    the jadaCorrectionSolver object
//! A_op            matrix A passed to jadaOp_create
//! B_op            matrix B passed to jadaOp_create
//! Qtil            projection vectors V passed to jadaOp_create
//! BQtil           projection vectors BV passed to jadaOp_create
//! sigma           (pos.!) shifts, -sigma[i], i in {1, ..., nvec} is passed to the jadaOp
//! res             JD residuals, e.g. rhs of the correction equations
//! resIndex        if not NULL, specifies permutation of the residual array to avoid unnecessary copying in the jada-algorithm
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! t               returns approximate solution vectors
//! iflag            a value > 0 indicates the number of systems that have not converged to the desired tolerance
void SUBR(jadaCorrectionSolver_run)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver,
                                    TYPE(const_op_ptr)    A_op,     TYPE(const_op_ptr)    B_op, 
                                    TYPE(const_mvec_ptr)  Qtil,     TYPE(const_mvec_ptr)  BQtil,
                                    const _ST_            sigma[],  TYPE(const_mvec_ptr)  res,      const int resIndex[], 
                                    const _MT_            tol[],    int                   maxIter,
                                    TYPE(mvec_ptr)        t,
                                    bool useIMGS,                   bool abortAfterFirstConvergedInBlock,
                                    int *                 iflag);
