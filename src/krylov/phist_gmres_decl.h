//! iteration status object for the pseudo-block GMRES
//! implemented in function gmres. To initialize the  
//! object before the first call, create the object   
//! using gmresState_create and adjust the parameters 
//! manually (the parameter maxBas must not be adjusted
//! after the call to _create as it determines the size
//! of the held data structures). Then call reset to   
//! set the initial vector. To restart from the        
//! previous iteration state, simply pass the object to
//! gmres again. To restart from an initial guess, use 
//! the reset function again. It is important to call  
//! reset before the first usage of the object in gmres
//! because the iteration will otherwise not start up  
//! correctly.
typedef struct TYPE(gmresState) {
  //! \name input and output args:
  //@{
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
  _MT_ tol; //! convergence tolerance for this system
  int  maxIters; //! maximum number of iterations allowed for this system
  int maxBas; //! maximum size of basis before restart
  int maxBasAllocated; //! maximum value that maxBas can take without deleting
                        //! and recreating the object.
  int ierr; //! error code returned by this GMRES instance
  //@}
  //! \name  internal data structures:
  //@{
  TYPE(mvec_ptr) V_; //! current basis
  TYPE(mvec_ptr) H_; //! current Hessenberg matrix from the Arnoldi process,
                     //! transformed to upper triangular form using Givens
                     //! rotations.
  _ST_ *cs_, *sn_;   //! cosine and sine terms for the Givens rotations
  _ST_ *rs_;
  
  int curDimV_; //! current size of the basis V
  int curIter_; //! number of iterations performed since last call to init
  //@}
} TYPE(gmresState);

typedef TYPE(gmresState)* TYPE(gmresState_ptr);

//! a simple GMRES implementation that works on several vectors simultaneously,
//! building a separate Krylov subspace for each of them. The iteration status 
//! is stored in a struct so that the process can be continued for some of the 
//! systems if one or more converge. This is not a block GMRES but a 'pseudo-  
//! block GMRES' as the former would build a single subspace for all rhs. It is
//! therefore OK to have the operator perform a different task for each vector 
//! column it is applied to, like in block JaDa: OP(X_j) = (A-s_jI)(I-VV').    
//! The computational steering is done by initializing the parameters in the   
//! gmresState structs. Individual ierr flags are contained within the structs,
//! and a global error code is put into the last arg ierr, as usual.
void SUBR(gmres)(TYPE(const_op_ptr) Op,
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B,
        TYPE(gmresState)* array_of_states,
        int* ierr);

//! delete gmresState object
void SUBR(gmresState_delete)(TYPE(mvec_ptr) S, int* ierr);

//! create a gmresState object for a given matrix/vector pair. The input parameters
//! will be set to default values and can be adjusted before calling gmres. 
//! The maxBas parameter cannot be adjusted afterwards as it determines the amount 
//! of memory allocated in this function. The map argument indicates the data dist-
//! ribution of vectors so we can create the basis vectors a priori.
SUBR(gmresState_create)(TYPE(gmresState_ptr)* state, const_map_ptr map, 
        int maxBas,int* ierr);

//! this function can be used to force a clean restart of the associated GMRES
//! solver. It is necessary to call this function before the first call to
//! gmres. The input starting vector x0 may be NULL, in that case this function
//! will generate a random initial guess. x0 does not have to be normalized in advance.
SUBR(gmresState_reset)(TYPE(gmresState_ptr) S, TYPE(const_mvec_ptr) x0,int *ierr);