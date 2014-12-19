
//! The feastCorrectionSolver uses carp_cg to calculate approximate solutions to a set 
//! of FEAST correction equations. It provides a simple interface and takes care of 
//! optimally pipelining the operations.
//! In each call to run(), numShifts linear systems with numBlocks right-hand sides
//! are solved simultaneously. Whenever a shift changes or the number of shifts/systems
//! changes, the object must be destroyed and rebuilt (This is just a first guess at
//! what may be a useful interface)
//! The systems to be solved are (sigma[j]I-A)(Xr[k]+i*Xi[k])=B[k],
//! where Xr, Xi and B are mvecs with blockSize columns, and sigma[j] are numShifts complex 
//! shifts.
//!
//! \todo for GHOST we should calculate in complex numbers directly but pass in a real-valued matrix (with complex shift)
typedef struct TYPE(feastCorrectionSolver)
{
  TYPE(const_sparseMat_ptr) A_;
  int numShifts_;       // number of shifts for which the system has to be solved
  _MT_ *sigma_r_;
  _MT_ *sigma_i_;
  int blockSize_;       // number of vectors in each system (#cols of rhs)
  TYPE(const_mvec_ptr) rhs_; // common right-hand side for all shifts sigma[j]
  linSolv_t method_;
  TYPE(carp_cgState_ptr) *carp_cgStates_;
  
} TYPE(feastCorrectionSolver);

typedef TYPE(feastCorrectionSolver)* TYPE(feastCorrectionSolver_ptr);

typedef TYPE(feastCorrectionSolver) const * TYPE(const_feastCorrectionSolver_ptr);


//! create a feastCorrectionSolver object. You have to specify the number of rhs 
//! per shift to be treated simultaneously (blockSize), the number of shifts 
//! (numShifts), and the complex shifts
void SUBR(feastCorrectionSolver_create)(TYPE(feastCorrectionSolver_ptr) *fCorrSolver, 
        TYPE(const_sparseMat_ptr) A, linSolv_t method,
        int blockSize, int numShifts,
        _MT_ shifts_r[], _MT_ shifts_i[],
        int *iflag);

//! delete a feastCorrectionSolver object
void SUBR(feastCorrectionSolver_delete)(TYPE(feastCorrectionSolver_ptr) fCorrSolver, int *iflag);


//! calculate approximate solutions to given set of FEAST correction equations
//!
//! arguments:
//! fCorrSolver    the feastCorrectionSolver object
//! rhs            rhs of the correction equations, rhs is the same for all
//!                systems corresponding to to shift[i] in create().
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! sol             returns approximate solution vectors, sol[i] belongs to shift[i] and rhs
//! iflag            a value > 0 indicates the number of systems that have not converged to the 
//!                 desired tolerance
//!
//! rhs and sol_r/sol_i must have the same number of columns (vectors), but it does *not* 
//! need to be identical to the blockSize specified in the create() function.
void SUBR(feastCorrectionSolver_run)(TYPE(feastCorrectionSolver_ptr) fCorrSolver,
                                    TYPE(const_mvec_ptr) rhs,
                                    const _MT_ tol, int maxIter,
                                    TYPE(mvec_ptr) sol_r[],
                                    TYPE(mvec_ptr) sol_i[],
                                    int *iflag);
