//! \addtogroup linear_solvers
//@{

//! \defgroup carp_cg CARP-CG row projection method for general (shifted) linear systems
//@{

// Pipelined CARP-CG solver for general shifted matrices sigma[j]I-A.
// The algorithm is CG on the minimum norm problem AA'x=b with SSOR pre-
// conditioning, implemented following the work of Bjoerck and Elfving (1979).
// The parallelization of the Kaczmarz sweeps is left to the kernel library
// (functions carp_setup, carp_sweep). We allow the special situation of a real matrix
// with complex shifts sigma = sigma_r + i*sigma_i, a real RHS vector and complex
// result x_r + i*x_i. Each state object may iterate on multiple systems with the
// same shift but different RHS (nvec>1)
//

//! state object for CARP-CG

//! the usage of carp_cg is very similar to that of blockedGMRES, except that we don't
//! need the complicated queuing of vectors in blockedGMRES.
typedef struct TYPE(carp_cgState) {

  //! \name output args:
  //@{
  int id; //! can be used to identify the system solved here (the shift to which this 
          //! iteration status belongs)
          //! This id is currently not used inside the carp_cg solver but might be used
          //! by the caller for reordering state objects etc. Initially, we just set
          //! id to i for object i created in _create().
  
  int numIter; //! number of iterations performed since last call to reset()
  _MT_ *normR;  //! stores current (implicit) residual norms

  int ierr; //! error code returned by this CARP-CG instance

  int *conv; //! set to 1 if system j has converged to the required tolerance

  //@}
  
  //! \name input data set by constructor
  //@{
  _MT_ sigma_r_; //! we're solving (sigma*I-A)x=b in this state object,
  _MT_ sigma_i_; //! with sigma = sigma_r + i*sigma_i
  int nvec_; //! number of RHS vectors for this shift

  TYPE(const_crsMat_ptr) A_;

  //@}
  //! \name set by reset() function
  //@{
  // rhs vector
  TYPE(const_mvec_ptr) b_;
  //@}
  //! \name internal CARP data structures
  //@{
  //! array with the 2-norm of each row of A-sigma*I (squared inverse, actually)
  _MT_* nrms_ai2i_;
  void* aux_; // work arg to carp_sweep (dep. on kernel lib)
  _MT_ omega_;// relaxation parameter
  //@}

  //! \name  internal CG data structures
  //@{
  TYPE(mvec_ptr) q_, r_, p_, z_; //! CG helper vectors, one column per RHS
  // imaginary parts, used if real A with complex shift
  TYPE(mvec_ptr) qi_, ri_, pi_, zi_;

  // scalars forming the Lanczos matrix, one for each RHS
  _ST_ *alpha_;
  _MT_ *alpha_i_, *beta_;

  _MT_ *normB_ ; //! two-norm of rhs vector (for stopping criteria)
  _MT_ *normR0_; //! stores initial (explicit) residual norms (for stopping criteria)
  _MT_ *normR_old;
  
  //@}
} TYPE(carp_cgState);

typedef TYPE(carp_cgState)* TYPE(carp_cgState_ptr);
typedef TYPE(carp_cgState) const * TYPE(const_cgState_ptr);

// constructor

//! Create an array of CG state objects to solve a set of numSys linear systems
//! (sigma[j]I-A)X=B. The systems all have the same rhs B, consisting of nvec  
//! columns.
void SUBR(carp_cgStates_create)(TYPE(carp_cgState_ptr) S_array[], 
        int nshifts, _MT_ sigma_r[], _MT_ sigma_i[],
        TYPE(const_crsMat_ptr) A, int numRhs, int* ierr);

//! destructor
void SUBR(carp_cgStates_delete)(TYPE(carp_cgState_ptr) S_array[], int numSys, int* ierr);

//! reset solver state

//! this function is used if the RHS vector changes but the matrix and shift
//! sigma stays the same. The solver is reset to start solving the new system
//! (sigma*I-A)x=rhs. If normsB==NULL, the two-norm of B is computed in S->normB_,
//! otherwise it is copied from the given pointer (length num_vectors of rhs).
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S,
                TYPE(const_mvec_ptr) rhs,
                _MT_* normsB,
                int *ierr);

//! CARP-CG iterations on all linear systems in the array.

//! To set the RHS vector, use reset() beforehand. If X is complex
//! or A and the shift sigma are real, X_i may be NULL. This function
//! performs up to maxIter CARP-CG steps on each of the numSys linear
//! systems and returns if all of them have either converged or failed.
//! Converged systems are 'locked', their x vectors are no longer updated
//! but for simplicity we keep doing operations like matvecs on them. It
//! is therefore advisable to not use too many columns per state object.
void SUBR(carp_cgStates_iterate)(
        TYPE(carp_cgState_ptr) S_array[], int numSys,
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ tol, int maxIter,
        int* ierr);

//@}

//@}
