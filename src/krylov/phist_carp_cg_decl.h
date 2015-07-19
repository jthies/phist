//! \addtogroup linear_solvers
//@{

//! \defgroup carp_cg CARP-CG row projection method for general (shifted) linear systems

//! Blocked CARP-CG solver for general shifted matrices A-sigma[j]I.
//! The algorithm is CGMN (CG on the minimum norm problem AA'x=b with SSOR pre-
//! conditioning, implemented following the work of Bjoerck and Elfving (1979)).
//! The parallelization of the Kaczmarz sweeps is left to the kernel library
//! (functions carp_setup, carp_sweep). 
//!
//! In order to be useful for the FEAST eigensolver, we allow the special situation of a real matrix
//! with complex shifts sigma = sigma_r + i*sigma_i, a real RHS vector and complex
//! result x_r + i*x_i. A state object may iterate on multiple systems with the
//! same or different shifts or RHS (nvec>1). 
//!
//! In order to be useful in the context of Jacobi-Davidson, we support adding projection vectors Q,
//! in which case the system solved is [A-sigma[j]I     V; trans(V)     0][x s]=[b 0]
//! The implementation is a bit simpler than for blocked GMRES because we can use a single state 
//! object for multiple shifts (no restart mechanism is needed for CG).
//@{

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

  int iflag; //! error code returned by this CARP-CG instance

  int *conv; //! set to 1 if system j has converged to the required tolerance

  //@}
  
  //! \name input data set by constructor
  //@{
  _MT_ *sigma_r_; //! we're solving (sigma*I-A)x=b in this state object,
  _MT_ *sigma_i_; //! with sigma = sigma_r + i*sigma_i
  int nvec_; //! number of RHS vectors for this shift
  bool rc_variant_;

  TYPE(const_sparseMat_ptr) A_;

  //@}
  //! \name set by reset() function
  //@{
  // rhs vector
  TYPE(const_mvec_ptr) b_;
  //@}
  //! \name internal CARP data structures
  //@{
  void* aux_; // work arg to carp_sweep (dep. on kernel lib)
  _MT_ *omega_;// relaxation parameter
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

//! Create a CG state object to solve a set of numSys linear systems
//! (A-sigma[j]I)X=B. For each column in B (=numSys), a shift must be
//! provided. sigma_i may be NULL if all shifts are real.
void SUBR(carp_cgState_create)(TYPE(carp_cgState_ptr) *S, 
        TYPE(const_sparseMat_ptr) A, 
        int numSys, _MT_ sigma_r[], _MT_ sigma_i[],
        int* iflag);

//! destructor
void SUBR(carp_cgState_delete)(TYPE(carp_cgState_ptr) S, int* iflag);

//! reset solver state

//! this function is used if the RHS vector changes but the matrix and shift
//! sigma stays the same. The solver is reset to start solving the new system
//! (sigma*I-A)x=rhs. If normsB==NULL, the two-norm of B is computed in S->normB_,
//! otherwise it is copied from the given pointer (length num_vectors of rhs).
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S,
                TYPE(const_mvec_ptr) rhs,
                _MT_* normsB,
                int *iflag);

//! CARP-CG iterations on all linear systems in the state object

//! To set the RHS vector, use reset() beforehand. If X is complex
//! or A and the shift sigma are real, X_i may be NULL. This function
//! performs up to maxIter CARP-CG steps on each of the linear systems
//! defined by the columns of x and b and the shifts passed to _create,
//! and returns if all of them have either converged or failed. If the
//! flag abortIfOneConverges is true, the function returns as soon as, 
//! well, guess what. Otherwise,
//! converged systems are 'locked', their x vectors are no longer updated
//! but for simplicity we keep doing operations like matvecs on them. It
//! is therefore advisable to not use too many columns per state object.
void SUBR(carp_cgState_iterate)(
        TYPE(carp_cgState_ptr) S,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        _MT_ tol, int maxIter, bool abortIfOneConverges,
        int* iflag);

//@}

//@}
