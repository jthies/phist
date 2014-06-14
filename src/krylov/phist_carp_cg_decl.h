// Pipelined CARP-CG solver for general shifted matrices sigma[j]I-A.
// The algorithm is CG on the minimum norm problem AA'x=b with SSOR pre-
// conditioning, implemented following the work of Bjoerck and Elfving (1979).
// The parallelization of the Kaczmarz sweeps is left to the kernel library
// (functions carp_setup, dkswp). We allow the special situation of a real matrix
// with complex shifts sigma = sigma_r + i*sigma_i, a real RHS vector and complex
// result x_r + i*x_i.
//

// the usage of carp_cg is very similar to that of pgmres, except that we don't
// need the complicated queuing of vectors in pgmres.
typedef struct TYPE(carp_cgState) {

//TODO - imaginary parts of vectors where needed.

  //! \name input and output args:
  //@{
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
          //! This id is currently not used inside the code, i.e. we assume state[i] belongs
          //! to X(:,i) everywhere. But it could be useful when reordering the state array 
          //! or printing debug info in a function which gets just one state object.
  _MT_ tol; //! convergence tolerance for this system
  
  _MT_ sigma_r; //! we're solving (sigma*I-A)x=b in this state object,
  _MT_ sigma_i; //! with sigma = sigma_r + i*sigma_i
  int nvec; //! number of RHS vectors for this shift

  int numIters; //! number of iterations performed since last call to reset()

  int ierr; //! error code returned by this CARP-CG instance

  //@}
  //! \name CARP data structures
  //@{
  //! array with the 2-norm of each row of A-sigma*I (squared inverse, actually)
  _MT_* norms_ai2i_;

  //! \name  internal CG data structures
  //@{
  TYPE(mvec_ptr) q_, r_, p_; //! CG helper vectors, one column per RHS
  TYPE(mvec_ptr) x0_; //! starting vector to compute the first residual vector r0
#ifndef IS_COMPLEX
  // imaginary parts, used if real A with complex shift
  TYPE(mvec_ptr) qi_, ri_, pi_;
  TYPE(mvec_ptr) x0i_; 
#endif  
  TYPE(mvec_ptr) b_; //! rhs to compute the first residual vector r0

  // scalars forming the Lanczos matrix, one for each RHS
  _MT_ *alpha_r_, *alpha_i_;
  _MT_ *beta_; 
  
  _MT_ *normR0_; //! stores initial (explicit) residual norms
  _MT_ *normR_; //! stores current (implicit) residual norms
  
  //@}
} TYPE(carp_cgState);

typedef TYPE(carp_cgState)* TYPE(carp_cgState_ptr);


typedef TYPE(carp_cgState) const * TYPE(const_cgState_ptr);

//! CARP-CG iterations on all linear systems simultaneously.
//! To set the RHS vector, use reset() beforehand. If A is complex
//! or the shift sigma is real, X_i may be NULL.
void SUBR(carp_cgStates_iterate)(TYPE(const_crsMat_ptr) A,
        TYPE(const_mvec_ptr) rhs, 
        TYPE(carp_cgState_ptr) S_array[], int numSys,
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        int* nIter, int* ierr);

//!
void SUBR(carp_cgStates_create)(TYPE(carp_cgState_ptr) S_array[], int numSys,
        const_map_ptr_t map, int numRhs, int* ierr);

//!
void SUBR(carp_cgStates_delete)(TYPE(carp_cgState_ptr) S_array[], int numSys, int* ierr);

//! this function can be used to force a clean restart of the associated CARP-CG
//! solver. It is necessary to call this function before the first call to
//! iterate(). The input starting vector x0 may be NULL, in that case this function
//! will generate a zero initial guess. x0 does not have to be normalized in advance.
//! x0i can be used to specify a complex initial vector (for real matrices with complex
//! shift). If A (and thus x0) is complex, x0i is ignored.
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S,
                TYPE(const_mvec_ptr) x0, TYPE(const_mvec_ptr) x0i, int *ierr);

