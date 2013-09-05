//! TODO doc me!
typedef struct TYPE(pfgmres_args)
  {

  // input/output

  TYPE(mvec_ptr) lhs;

  // input arguments

  TYPE(const_mvec_ptr) rhs;

  // pseudo-thread id (for instance vector column on which to work)
  int t_id;
  // the first three arguments allow
  // the pfgmres function to perform
  // (A-sigma*I)X and M\X, we just  
  // need to initialize the task    
  // buffer correctly.
  taskBuf_t *taskBuf;
  
  int op_AX; // flag to request A*x from the task buffer
  int op_RPRECX; // flag to request Mr\x from the task buffer
  int op_nrm2;   // flag to request ||v||_2
  int op_vTw;    // flag to request C=V'W
  
  _ST_ shift; // the operator used in the algorithm is A + shift*I
  
  // GMRES parameters
  int max_iters;
  int num_blocks;
  _MT_ tol;

  // output arguments
  int num_iters;
  int ierr;
  
  } TYPE(pfgmres_args);

//! Pipelined Generalized Minimum Residual method.
//! this variant of GMRES is intended for the parallel
//! solution of multiple linear systems of the form
//! (A-sigma_j I)x_j=b_j, j=1..num_sys. However, any 
//! kind of operator can be used. The number of
//! columns in X and B must be num_sys (This is not a
//! block variant of GMRES).
void SUBR(pfgmres)(TYPE(pfgmres_args)* args);
