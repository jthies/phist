// Pipelined CARP-CG solver for general shifted matrices A-sigma[j]B.
// The algorithm is CG on the minimum norm problem AA'x=b with SSOR pre-
// conditioning, implemented following the work of Bjoerck and Elfving (1979).
// The parallelization of the Kaczmarz sweeps is left to the kernel library
// (function carp_fb).
//

// the usage of carp_cg is very similar to that of pgmres, except that we don't
// need the complicated queuing of vectors in pgmres.
typedef struct TYPE(carp_cgState) {
  //! \name input and output args:
  //@{
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
          //! This id is currently not used inside the code, i.e. we assume state[i] belongs
          //! to X(:,i) everywhere. But it could be useful when reordering the state array 
          //! or printing debug info in a function which gets just one state object.
  _MT_ tol; //! convergence tolerance for this system
  
  _ST_ sigma; //! we're solving (A-sigma*I)x=b in this state object

  int ierr; //! error code returned by this CARP-CG instance

  int totalIter; //! counts number of iterations (also over restarts)
  //@}
  //! \name CARP data structures
  //@{
  //! vector with the 2-norm of each row of A-sigma*I (squared)
  TYPE(mvec_ptr) norms_ai2_;

  //! \name  internal CG data structures
  //@{
  TYPE(mvec_ptr) q_, r_, p_; //! this instance operates on column 'id' of these vectors only.
  TYPE(mvec_ptr) x0_; //! starting vector to compute the first residual vector r0
  TYPE(mvec_ptr) b_; //! rhs to compute the first residual vector r0
  _MT_ alpha_, beta_; // scalars forming the Lanczos (cf. numerics textbook for CG algorithm)
  
  _MT_ normR0_; //! stores initial (explicit) residual norm
  _MT_ normR_; //! stores current (implicit) residual norm
  
  int maxIters_; //! maximum number of iterations allowed

  //@}
} TYPE(carp_cgState);

typedef TYPE(carp_cgState)* TYPE(carp_cgState_ptr);

typedef TYPE(carp_cgState) const * TYPE(const_cgState_ptr);

//!
void SUBR(carp_cgStates_iterate)(TYPE(const_crsMat_ptr) A,
        TYPE(carp_cgState_ptr) S_array[], TYPE(mvec_ptr) X,
        int* nIter, int* ierr);

//!
void SUBR(carp_cgStates_create)(TYPE(carp_cgState_ptr) S_array[], int numSys,
        const_map_ptr_t map, int maxIters, int* ierr);

//!
void SUBR(carp_cgStates_delete)(TYPE(carp_cgState_ptr) S_array[], int numSys, int* ierr);

//! this function can be used to force a clean restart of the associated CARP-CG
//! solver. It is necessary to call this function before the first call to
//! pcg. The input starting vector x0 may be NULL, in that case this function
//! will generate a random initial guess. x0 does not have to be normalized in advance.
//! The input crsMat A and RHS B may also be NULL, meaning 'keep old one', but not on 
//! the first call to reset. If one of the RHS vectors changes between calls to pcg, 
//! reset with the new rhs should be called for that cgState, otherwise a messed up 
//! Krylov sequence will result and the convergence criterion will not be consistent. If
//! the matrix A changes, *all* states must be reset. If some the shift in JaDa changes,
//! reset should be called with A=NULL (A does not change) for that particular system.
void SUBR(carp_cgState_reset)(TYPE(carp_cgState_ptr) S,
                TYPE(const_crsMat_ptr) A, TYPE(const_mvec_ptr) b,
                TYPE(const_mvec_ptr) x0, int *ierr);
