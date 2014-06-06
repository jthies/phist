
#ifdef IS_DOUBLE
typedef Zmvec_t Dfeast_mvec_t;
typedef d_complex_t Dfeast_shift_t;
#else
typedef Cmvec_t Sfeast_mvec_t;
typedef s_complex_t Sfeast_shift_t;
typedef 
#endif

//!
//! The feastCorrectionSolver uses carp_cg to calculate approximate solutions to a set 
//! of FEAST correction equations. It provides a simple interface and takes care of 
//! optimally pipelining the operations.
//! In each call to solve(), numShifts linear systems with numBlocks right-hand sides
//! are solved simultaneously. Whenever a shift changes or the number of shifts/systems
//! changes, the object must be destroyed and rebuilt (This is just a first guess at
//! what may be a useful interface)
typedef struct TYPE(feastCorrectionSolver)
{
// TODO
} TYPE(feastCorrectionSolver);

typedef TYPE(feastCorrectionSolver)* TYPE(feastCorrectionSolver_ptr);

typedef TYPE(feastCorrectionSolver) const * TYPE(const_feastCorrectionSolver_ptr);


//! create a feastCorrectionSolver object. You have to specify the number of rhs 
//! per shift (blockSize), the number of shifts (numShifts), and the complex shifts

void SUBR(feastCorrectionSolver_create)(TYPE(feastCorrectionSolver_ptr) *fCorrSolver, 
        TYPE(const_crsMat_ptr) A, 
        int blockSize, int numShifts,
        TYPE(feast_shift) shifts[],
        int *ierr);

//! delete a feastCorrectionSolver object
void SUBR(feastCorrectionSolver_delete)(TYPE(feastCorrectionSolver_ptr) fCorrSolver, int *ierr);


//! calculate approximate solutions to given set of FEAST correction equations
//!
//! arguments:
//! fCorrSolver    the feastCorrectionSolver object
//! rhs            rhs of the correction equations, rhs is the same for all
//!                to shift[i] in create()
//! tol             desired accuracy (gmres residual tolerance) of the individual systems
//! maxIter         maximal number of iterations after which individial systems should be aborted
//! sol             returns approximate solution vectors, sol[i] belongs to shift[i] and rhs
//! ierr            a value > 0 indicates the number of systems that have not converged to the 
//!                 desired tolerance
void SUBR(feastCorrectionSolver_run)(TYPE(feastCorrectionSolver_ptr) fCorrSolver,
                                    TYPE(feast_mvec) const * const rhs,
                                    const _MT_ tol, int maxIter,
                                    TYPE(feast_mvec) * sol[],
                                    int *ierr);
