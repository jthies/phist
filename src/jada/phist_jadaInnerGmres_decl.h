//! TODO: update description, partially outdated, refers to phist_gmres rather thand this variant

//! iteration status object for the pseudo-block GMRES
//! iteration we currently use as approximate solver  
//! for the correction equation. To initialize the    
//! object before the first call, create the object   
//! using gmresState_create and adjust the parameters 
//! manually (the parameter maxBas must not be adjusted
//! after the call to _create as it determines the size
//! of the held data structures). Then call reset to   
//! set the rhs and initial vector. To restart from the
//! previous iteration state, simply pass the object to
//! gmres again. To restart from an initial guess, use 
//! the reset function again. It is important to call  
//! reset before the first usage of the object in gmres
//! because the iteration will otherwise not start up  
//! correctly. See test/krylov/TestGmres for examples  
//! of using this object to build up a restarted GMRES.
//!                                                    
//! In order to utilize 'scattered views' when e.g.    
//! applying A*X to vectors from different gmresStates,
//! we need to draw all the bases V from a single block
//! of memory, Vglob. The way we implement this is to  
//! always create and delete an entire array of states.
//! not all of them have to be passed to functions like
//! iterate or updateSol, but they must be created and 
//! deleted together so that the memory handling is    
//! clear.                                             
typedef struct TYPE(jadaInnerGmresState) {
  //! \name input and output args:
  //@{
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
          //! This id is currently not used inside the code, i.e. we assume state[i] belongs
          //! to X(:,i) everywhere. But it could be useful when reordering the state array 
          //! or printing debug info in a function which gets just one state object.
  _MT_ tol; //! convergence tolerance for this system

  int ierr; //! error code returned by this GMRES instance

  int totalIter; //! counts number of iterations (also over restarts)
  //@}
  //! \name  internal data structures
  //@{
  TYPE(mvec_ptr) *V_; //! array of mvecs because row-wise storage is not well suited here
  TYPE(mvec_ptr) *AV_; //! ""
  void *glob_unused_mvecs_;
  void *glob_used_mvecs_;
  TYPE(sdMat_ptr) H_; //! memory block in which Hessenberg matrix from the Arnoldi process
                     //! is built up,, transformed to upper triangular form using 
                     //! Givens rotations.
  TYPE(mvec_ptr) x0_; //! starting vector to compute the first residual vector r0
  TYPE(mvec_ptr) b_; //! rhs to compute the first residual vector r0
  _MT_ *cs_;         //! terms c of the Givens rotation
  _ST_ *sn_;         //! terms s of the Givens rotation
  _ST_ *rs_;
  
  _MT_ normR0_; //! stores initial (explicit) residual norm
  _MT_ normR_; //! stores current (implicit) residual norm
  
  int curDimV_; //! current size of the basis V

  int maxBas_; //! maximum size of basis before restart

  //@}
} TYPE(jadaInnerGmresState);

typedef TYPE(jadaInnerGmresState)* TYPE(jadaInnerGmresState_ptr);

typedef TYPE(jadaInnerGmresState) const * TYPE(const_jadaInnerGmresState_ptr);

//!                                                                                     
//! a simple GMRES implementation that works on several vectors simultaneously,         
//! building a separate Krylov subspace for each of them. The iteration status          
//! is stored in a struct so that the process can be continued for some of the          
//! systems if one or more converge. This is not a block GMRES but a 'pseudo-           
//! block GMRES' as the former would build a single subspace for all rhs. It is         
//! therefore OK to have the operator perform a different task for each vector          
//! column it is applied to, like in block JaDa: OP(X_j) = (A-s_jI)(I-VV').             
//! (Note: a real block variant should be possible, since A-s_jI gives the *SAME* Krylov-subspace as A)
//! The computational steering is done by initializing the parameters in the            
//! gmresState structs. The iteration stops as soon as one system converges or          
//! reaches maxBas iterations (NOTE: no restarting is implemented, so the maximum       
//! number of iterations is determined by the number of basis vectors allocated).       
//! The user can then continue running GMRES on the others after appropriately          
//! reordering the rhs vectors and array of states.                                     
//!                                                                                     
//! This function does *not* compute the solution to the original system AX=B for you.  
//! For any of the systems (wether converged or not) the function gmresState_updateSol  
//! can be used to compute the solution from the state object and the approximation     
//! passed to the previous reset() call.                                                
//!                                                                                     
//! Individual ierr flags are contained within the structs,                             
//! and a global error code is put into the last arg ierr, as usual.                    
//!                                                                                     
//! if system j has converged, S_array[j]->ierr will be set to 0 on output.             
//! Otherwise it will be set to 1 (not yet converged) or 2 (max iters exceeded),        
//! or a negative value if an error occurred related to this particular system.         
//! The global ierr flag will then be set to -1. (0 for "someone converged" and +1 for  
//! someone reached max iters")                                                         
//!                                                                                     
void SUBR(jadaInnerGmresStates_iterate)(TYPE(const_op_ptr) Op,
        TYPE(jadaInnerGmresState_ptr) S_array[], int numSys,
        int* nIter, int* ierr);

//! create an array of gmresState objects. The method's input parameters
//! will be set to default values and can be adjusted before calling reset/iterate. 
//! The maxBas parameter cannot be adjusted afterwards as it determines the amount 
//! of memory allocated in this function. The map argument indicates the data dist-
//! ribution of vectors so we can create the basis vectors a priori.
//! The array of pointers must be allocated beforehand, but the individual structs 
//! are allocated by this method.
void SUBR(jadaInnerGmresStates_create)(TYPE(jadaInnerGmresState_ptr) S_array[], int numSys,
        const_map_ptr_t map, int maxBas, int* ierr);

//! delete an set of gmresState objects. Only the individual structs are destroyed,
//! The csller has to delete the array and nullify it.
void SUBR(jadaInnerGmresStates_delete)(TYPE(jadaInnerGmresState_ptr) S_array[], int numSys, int* ierr);

//! this function can be used to force a clean restart of the associated GMRES
//! solver. It is necessary to call this function before the first call to
//! gmres. The input starting vector x0 may be NULL, in that case this function
//! will generate a random initial guess. x0 does not have to be normalized in advance.
//! The input RHS may also be NULL, meaning 'keep old RHS', but not on the first call to
//! reset. If one of the RHS vectors changes between calls to gmres, reset with the new
//! rhs should be called for that gmresState, otherwise a messed up Krylov sequence will
//! result and the convergence criterion will not be consistent.
void SUBR(jadaInnerGmresState_reset)(TYPE(jadaInnerGmresState_ptr) S,
                TYPE(const_mvec_ptr) b,
                TYPE(const_mvec_ptr) x0, int *ierr);

//! For each of the state objects i passed in, update the current approximation x(:,i) using 
//! the basis V and the projection coefficients. This should be done for a system that 
//! indicates convergence after iterate, but it can also be done to get an intermediate 
//! solution. The function is 'vectorized' in the same way as iterate, so an array of states 
//! and multivector x can be passed in.
//! If not NULL, also outputs A*x through projections coefficients and A*V
void SUBR(jadaInnerGmresStates_updateSol)(TYPE(jadaInnerGmresState_ptr) S_array[], int numSys,
                TYPE(mvec_ptr) x, TYPE(mvec_ptr) Ax, _MT_ *resNorm, bool scaleSolutionToOne, int* ierr);
