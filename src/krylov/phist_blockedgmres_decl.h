//! \defgroup linear_solvers Iterative methods for linear systems
//@{

//! \defgroup blockedGMRES blocked GMRES solver for general linear systems
//@{

//! gmres state object

//! iteration status object for the pipelined GMRES
//! iteration we currently use as approximate solver
//! for the correction equation. To initialize the
//! object before the first call, create the object
//! using gmresState_create and adjust the parameters 
//! manually (the parameter maxDimV must not be adjusted
//! after the call to _create as it determines the size
//! of the held data structures). Then call reset to
//! set the rhs and initial vector. To restart from the
//! previous iteration state, simply pass the object to
//! gmres again. To restart from an initial guess, use 
//! the reset function again. It is important to call
//! reset before the first usage of the object in gmres
//! because the iteration will otherwise not start up
//! correctly. See test/krylov/TestBlockedGMRES for examples
//! of using this object to build up a restarted GMRES.
typedef struct TYPE(blockedGMRESstate)
{
  //! \name input and output args:
  //@{
  int id;               //! used to identify the system solved, set in blockedGMRESSates_create, don't modify!
  _MT_ tol;             //! convergence tolerance for this system (can be adjusted any time)
  int status;           //! -2: not initialized, -1: resetted, 0: converged, 1: not yet converged, 2: max iters exceeded
  int totalIter;        //! counts the total number of iterations (also over restarts)
  //@}
  //! \name  internal data structures
  //@{
  int lastVind_;        //! index of the last element in the V-buffer
  int curDimV_;         //! current dimension of the subspace in V_
  TYPE(mvec_ptr) b_;    //! right hand side vector of the system Ax=b
  TYPE(sdMat_ptr) H_;   //! Hessenberg-matrix from the Arnoldi-process rotated to upper triangular form
  _ST_ *cs_;            //! terms c of the Givens rotations (c s; -s' c)
  _ST_ *sn_;            //! terms s of the Givens rotations (c s; -s' c)
  _ST_ *rs_;            //! rotated projected residual (e.g. ...*q_3*q_2*q_1*e1 )
  _MT_ normR0_;         //! initial (explicit) residual norm
  _MT_ normR_;          //! current (implicit) residual norm
  _ST_ prevBeta_;       //! previous secondary diagonal entry for the MINRES variant (from the unrotated Hessenberg (here tridiagonal) matrix)

  void *Vbuff;          //! ring buffer for the subspaces V
  //@}
} TYPE(blockedGMRESstate);

typedef TYPE(blockedGMRESstate)* TYPE(blockedGMRESstate_ptr);

typedef TYPE(blockedGMRESstate) const * TYPE(const_blockedGMRESstate_ptr);

//!
//! a simple GMRES implementation that works on several vectors simultaneously,
//! building a separate Krylov subspace for each of them. The iteration status
//! is stored in a struct so that the process can be continued for some of the
//! systems if one or more converge. This is not a block GMRES but a 'pseudo- 
//! block GMRES' as the former would build a single subspace for all rhs. It is
//! therefore OK to have the operator perform a different task for each vector
//! column it is applied to, like in block JaDa: OP(X_j) = (A-s_jI)(I-VV').
//! (Note: a real block variant should be possible, since A-s_jI gives the *SAME* Krylov-subspace as A,
//! but it would be difficult to add a new system to an existing (generalized) krylov-subspace)
//! The computational steering is done by initializing the parameters in the
//! gmresState structs. The iteration stops as soon as one system converges or
//! reaches maxDimV iterations (NOTE: no restarting is implemented, so the maximum
//! number of iterations is determined by the number of basis vectors allocated).
//! The user can then continue running GMRES on the others after appropriately
//! reordering the rhs vectors and array of states.
//!
//! This function does *not* compute the solution to the original system AX=B for you.
//! For any of the systems (wether converged or not) the function gmresState_updateSol
//! can be used to compute the solution from the state object and the approximation
//! passed to the previous reset() call.
//!
//! Individual status flags are contained within the structs,
//! and a global error code is put into the last arg ierr, as usual.
//!
//! if system j has converged, S_array[j]->status will be set to 0 on output.
//! Otherwise it will be set to 1 (not yet converged) or 2 (max iters exceeded),
//! or a negative value if an error occurred related to this particular system.
//! The global ierr flag will then be set to -1. (0 for "someone converged" and +1 for
//! someone reached max iters")
//! \warning you cannot mix together states from different calls to blockedGMRESstates_create!
//!
void SUBR( blockedGMRESstates_iterate ) (TYPE(const_op_ptr) Op, TYPE(blockedGMRESstate_ptr) S_array[], int numSys, int* nIter, bool useIMGS, int* ierr);

//!
//! create an array of gmresState objects. The method's input parameters
//! will be set to default values and can be adjusted before calling reset/iterate. 
//! The maxBas parameter cannot be adjusted afterwards as it determines the amount 
//! of memory allocated in this function. The map argument indicates the data dist-
//! ribution of vectors so we can create the basis vectors a priori.
//! The array of pointers must be allocated beforehand, but the individual structs 
//! are allocated by this method.
//!
void SUBR( blockedGMRESstates_create ) (TYPE(blockedGMRESstate_ptr) S_array[], int numSys, const_map_ptr_t map, int maxBas, int* ierr);

//!
//! delete an set of gmresState objects. Only the individual structs are destroyed,
//! The csller has to delete the array and nullify it.
//! \warning you cannot delete individual states, but must pass the whole array created with blockedGMRESstates_create
//!
void SUBR( blockedGMRESstates_delete ) (TYPE(blockedGMRESstate_ptr) S_array[], int numSys, int* ierr);

//!
//! this function can be used to force a clean restart of the associated GMRES
//! solver. It is necessary to call this function before the first call to
//! blockedGMRES_iterate. The input starting vector x0 may be NULL, in that case this function
//! will generate a random initial guess. x0 does not have to be normalized in advance.
//! The input RHS may also be NULL, meaning 'keep old RHS', but not on the first call to
//! reset. If one of the RHS vectors changes between calls to gmres, reset with the new
//! rhs should be called for that gmresState, otherwise a messed up Krylov sequence will
//! result and the convergence criterion will not be consistent.
//!
void SUBR( blockedGMRESstate_reset ) (TYPE(blockedGMRESstate_ptr) S, TYPE(const_mvec_ptr) b, TYPE(const_mvec_ptr) x0, int *ierr);

//!
//! For each of the state objects i passed in, update the current approximation x(:,i) using 
//! the basis V and the projection coefficients. This should be done for a system that 
//! indicates convergence after iterate, but it can also be done to get an intermediate 
//! solution. The function is 'vectorized' in the same way as iterate, so an array of states 
//! and multivector x can be passed in.
//!
void SUBR( blockedGMRESstates_updateSol ) (TYPE(blockedGMRESstate_ptr) S_array[], int numSys, TYPE(mvec_ptr) x, _MT_ *resNorm, bool scaleSolutionToOne, int* ierr);

//@}
//@}