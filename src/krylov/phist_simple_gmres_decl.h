//! iteration status object for the pseudo-block GMRES
//! implemented in function simple_gmres.
typedef struct TYPE(gmresStatus) {
  // input and output args:
  int id; //! can be used to identify the system solved here (the column to which this 
          //! iteration status belongs)
  _MT_ tol; //! convergence tolerance for this system
  int  max_iters; //! maximum number of iterations allowed for this system
  int ierr; //! error code returned by this GMRES instance

  // internal data structures:
  TYPE(mvec_ptr) V_; // current basis
  TYPE(mvec_ptr) H_; // current Hessenberg matrix from the Arnoldi process
  _ST_ *cs_, *sn_;   // cosine and sine terms for the Givens rotations
  _ST_ *rs_;

} TYPE(gmresStatus);

//! a simple GMRES implementation that works on several vectors simultaneously,
//! building a separate Krylov subspace for each of them. The iteration status 
//! is stored in a struct so that the process can be continued for some of the 
//! systems if one or more converge. This is not a block GMRES but a 'pseudo-  
//! block GMRES' as the former would build a single subspace for all rhs. It is
//! therefore OK to have the operator perform a different task for each vector 
//! column it is applied to, like in block JaDa: OP(X_j) = (A-s_jI)(I-VV').    
//! The computation steering is done by initializing the parameters in the     
//! gmres_t struct.
void SUBR(simple_gmres)(TYPE(const_op_ptr) Op,
        TYPE(mvec_ptr) X,
        TYPE(const_mvec_ptr) B,
        int* ierr);
