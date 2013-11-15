//! iteration status object for the pseudo-block GMRES
//! implemented in function simple_gmres. To initialize
//! the object before the first call, create the object
//! using gmresState_create and adjust the parameters 
//! manually (the parameter maxBas must not be adjusted
//! after the call to _create as it determines the size
//! of the held data structures).                      
typedef struct TYPE(gmresState) {
  //! \name input and output args:
  //@{
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
  _MT_ tol; //! convergence tolerance for this system
  int  max_iters; //! maximum number of iterations allowed for this system
  int maxBas; //! maximum size of basis before restart
  int ierr; //! error code returned by this GMRES instance
  //@}
  //! \name  internal data structures:
  //@{
  TYPE(mvec_ptr) V_; //! current basis
  TYPE(mvec_ptr) H_; //! current Hessenberg matrix from the Arnoldi process
  _ST_ *cs_, *sn_;   //! cosine and sine terms for the Givens rotations
  _ST_ *rs_;
  
  int curDimV_; // current size of the basis V
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
void SUBR(simple_gmres)(TYPE(const_op_ptr) Op,
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B,
        TYPE(gmresState)* array_of_states,
        int* ierr);

//! create a gmresState object for a given matrix/vector pair. The input parameters
//! will be set to default values and can be adjusted before calling simple_gmres. 
//! The maxBas parameter cannot be adjusted afterwards as it determines the amount 
//! of memory allocated in this function. The map argument indicates the data dist-
//! ribution of vectors so we can create the basis vectors a priori.
SUBR(gmresState_create)(TYPE(gmresState_ptr)* state, const_map_ptr map, 
        int maxBas,int* ierr);

//! this function can be used to force a clean restart of the associated GMRES
//! solver. It is not necessary to call this function before the first call to
//! simple_gmres, but it can be used in that case just to set the starting    
//! vector.
SUBR(gmresState_init)(TYPE(gmresState_ptr) S, TYPE(const_mvec_ptr) x0,int *ierr);
