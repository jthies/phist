
// these are some thoughts I wrote down some time ago, I keep
// them here for the moment to look at later

//! in Jacobi-Davidson, we want ot perform an iterative solve
//! (GMRES to start with) for the correction equation. The   
//! corresponding operator we put in the GMRES solver is     
//!                                                          
//! (I-VV') M\(A-sigma_j I),                                 
//!                                                          
//! where V = [Q v] is the basis                             
//! of converged eigenvectors combined with the current      
//! search space, and M is some preconditioner for a nearby  
//! matrix A-tau*I). We want the GMRES method to run in para-
//! lel for several shifts sigma_j, which we achieve by using
//! the task-buffer programming model (cf. examples/sched/...
//! main_task_model.c for a simple example). For now we let  
//! the jadaOp handle the orthogonalization wrt. V, but in   
//! the end we can hopefully use some bordered preconditio-  
//! ning (like in HYMLS) for the matrix                      
//!                                                          
//! A-tau*I  V                                               
//!     V'   0                                               
//!                                                          
//! There are in principle two options, use (A-sigma_jI) as  
//! operator and "(I-VV')M\" as left preconditioner, or use  
//! (I-VV')(A-sigma_jI) as operator and M\ as right precond. 
//! The first is more suitable for the case where the pre-   
//! conditioner can handle the border, the second if flexible
//! ible GMRES has to be used, e.g. in case of nested ite-   
//! rations on the separators in DSC. Combining the right    
//! preconditioner with the projection should also be fine,  
//! the operator is then (A-sigma_j I) M\(I-VV'), I think    
//! that should give the same results more or less.
//!                                                          


// private struct to keep all the pointers we need in order to
// apply the operator.
typedef struct TYPE(jadaOp_data)
  {
  TYPE(const_op_ptr) A_op;
  TYPE(const_op_ptr) B_op;
  TYPE(const_mvec_ptr) V;
  _ST_ sigma;
  } TYPE(jadaOp_data);

//
void SUBR(jadaOp_apply)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
                        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;// not implemented
  return;
  }


//
void SUBR(jadaOp_create)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) V, _ST_ shift, 
                         TYPE(op_ptr)* jdOp, int *ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  *jdOp=(TYPE(op_ptr))malloc(sizeof(TYPE(op_ptr)));
  TYPE(jadaOp_data) *myOp=(TYPE(jadaOp_data)*)malloc(sizeof(TYPE(jadaOp_data)));

  myOp->A_op=A_op;
  myOp->B_op=B_op;
  myOp->V=V;
  myOp->sigma=shift;//TODO - define sign of sigma/shift
  
  (*jdOp)->A=(void*)myOp;
  (*jdOp)->apply = &SUBR(jadaOp_apply);
  return;
  }

