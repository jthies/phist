//! create a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_create)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver, int pgmresBlockDim, const_map_ptr_t map, int pgmresMaxBase, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

}

//! delete a jadaCorrectionSolver object
void SUBR(jadaCorrectionSolver_delete)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver, int *ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

}


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
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! t               returns approximate solution vectors
//! ierr            a value > 0 indicates the number of systems that have not converged to the desired tolerance
void SUBR(jadaCorrectionSolver_run)(TYPE(jadaCorrectionSolver_ptr) jdCorrSolver,
                                    TYPE(const_op_ptr)    A_op,     TYPE(const_op_ptr)    B_op, 
                                    TYPE(const_mvec_ptr)  Qtil,     TYPE(const_mvec_ptr)  BQtil,
                                    const _ST_            sigma[],  TYPE(const_mvec_ptr)  res,
                                    const _MT_            tol[],    int                   maxIter,
                                    TYPE(mvec_ptr)        t,        int *                 ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

}
