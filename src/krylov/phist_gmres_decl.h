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
typedef struct TYPE(gmresState) {
  //! \name input and output args:
  //@{
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
  _MT_ tol; //! convergence tolerance for this system
  int  maxIters; //! maximum number of iterations allowed for this system
  int maxBas; //! maximum size of basis before restart

  int ierr; //! error code returned by this GMRES instance
  //@}
  //! \name  internal data structures
  //@{
  TYPE(const_mvec_ptr) B_; //! for which RHS is this GMRES being run?
  TYPE(mvec_ptr) X0_; //! initial vector passed in to most recent reset() call (will be 
                      //! copied to make sure the user doesn't change/delete it, we need 
                      //! it in updateSol)
  TYPE(mvec_ptr) V_; //! current basis
  TYPE(mvec_ptr) H_; //! current Hessenberg matrix from the Arnoldi process,
                     //! transformed to upper triangular form using Givens
                     //! rotations.
  _ST_ *cs_, *sn_;   //! cosine and sine terms for the Givens rotations
  _ST_ *rs_;
  
  _MT_ normB_; //! stores norm of RHS so it doesn't have to be recomputed when
               //! when continuing.

  _MT_ normR_; //! stores current (implicit) residual norm
  
  int maxBasAllocated_; //! maximum value that maxBas can take without deleting
                        //! and recreating the object.
  int curDimV_; //! current size of the basis V
  int curIter_; //! number of iterations performed since last call to reset with new RHS

  //@}
} TYPE(gmresState);

typedef TYPE(gmresState)* TYPE(gmresState_ptr);

typedef TYPE(gmresState) const * TYPE(const_gmresState_ptr);

//!                                                                                     
//! a simple GMRES implementation that works on several vectors simultaneously,         
//! building a separate Krylov subspace for each of them. The iteration status          
//! is stored in a struct so that the process can be continued for some of the          
//! systems if one or more converge. This is not a block GMRES but a 'pseudo-           
//! block GMRES' as the former would build a single subspace for all rhs. It is         
//! therefore OK to have the operator perform a different task for each vector          
//! column it is applied to, like in block JaDa: OP(X_j) = (A-s_jI)(I-VV').             
//! The computational steering is done by initializing the parameters in the            
//! gmresState structs. The iteration stops as soon as one system converges or          
//! reaches maxBas iterations (NOTE: no restarting is implemented, so the maximum       
//! number of iterations is determined by the number of basis vectors allocated).       
//! The user can then continue running GMRES on the others after appropriately          
//! reordering the rhs vectors and array of states.                                     
//!                                                                                     
//! This function does *not* compute the solution to the original system AX=B for you.  
//! For any of the systems (wether converged or not) the function gmresState_getSol     
//! can be used to compute the solution from the state object.                          
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
void SUBR(gmresState_iterate)(TYPE(const_op_ptr) Op,
        TYPE(gmresState_ptr) S_array[], int numSys,
        int* ierr);

//! create a gmresState object for a given matrix/vector pair. The input parameters
//! will be set to default values and can be adjusted before calling gmres. 
//! The maxBas parameter cannot be adjusted afterwards as it determines the amount 
//! of memory allocated in this function. The map argument indicates the data dist-
//! ribution of vectors so we can create the basis vectors a priori.
void SUBR(gmresState_create)(TYPE(gmresState_ptr) S_array[], const_map_ptr_t map,
        int maxBas, int* ierr);

//! delete gmresState object
void SUBR(gmresState_delete)(TYPE(gmresState_ptr) S, int* ierr);

//! this function can be used to force a clean restart of the associated GMRES
//! solver. It is necessary to call this function before the first call to
//! gmres. The input starting vector x0 may be NULL, in that case this function
//! will generate a random initial guess. x0 does not have to be normalized in advance.
//! The input RHS may also be NULL, meaning 'keep old RHS', but not on the first call to
//! reset. If one of the RHS vectors changes between calls to gmres, reset with the new
//! rhs should be called for that gmresState, otherwise a messed up Krylov sequence will
//! result and the convergence criterion will not be consistent.
void SUBR(gmresState_reset)(TYPE(gmresState_ptr) S,
                TYPE(const_mvec_ptr) b,
                TYPE(const_mvec_ptr) x0, int *ierr);

//! update the current approximation x using the basis V and the projection coefficients.  
//! This should be done for a system that indicates convergence after iterate, but it can  
//! also be done to get an intermediate solution. The function is 'vectorized' in the same 
//! way as iterate, so an array of states and multivector x can be passed in.
void SUBR(gmresState_updateSol)(TYPE(gmresState_ptr) S_array[], TYPE(mvec_ptr) x, int* ierr);
