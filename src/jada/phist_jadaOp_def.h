
// these are some thoughts I wrote down some time ago, I keep
// them here for the moment to look at later

//! in Jacobi-Davidson, we want ot perform an iterative solve
//! (GMRES to start with) for the correction equation. The   
//! corresponding operator we put in the GMRES solver is     
//!                                                          
//! (I-BVV') M\(A-sigma_j B),                                 
//!                                                          
//! where V = [Q v] is the basis                             
//! of converged eigenvectors combined with the current      
//! search space, and M is some preconditioner for a nearby  
//! matrix A-tau*B). We want the GMRES method to run in para-
//! lel for several shifts sigma_j, which we achieve by using
//! the task-buffer programming model (cf. examples/sched/...
//! main_task_model.c for a simple example). For now we let  
//! the jadaOp handle the orthogonalization wrt. V, but in   
//! the end we can hopefully use some bordered preconditio-  
//! ning (like in HYMLS) for the matrix                      
//!                                                          
//! A-tau*B  BV                                               
//!     V'   0                                               
//!                                                          
//! There are in principle two options, use (A-sigma_jI) as  
//! operator and "(I-BVV')M\" as left preconditioner, or use  
//! (I-BVV')(A-sigma_jB) as operator and M\ as right precond. 
//! The first is more suitable for the case where the pre-   
//! conditioner can handle the border, the second if flexible
//! ible GMRES has to be used, e.g. in case of nested ite-   
//! rations on the separators in DSC. Combining the right    
//! preconditioner with the projection should also be fine,  
//! the operator is then (A-sigma_j B) M\(I-VV'B), I think    
//! that should give the same results more or less.
//!                                                          


// private struct to keep all the pointers we need in order to apply the operator.
typedef struct TYPE(jadaOp_data)
{
  TYPE(const_op_ptr)    A_op;   // operator of the general matrix A
  TYPE(const_op_ptr)    B_op;   // operator of the hpd. matrix B, assumed I when NULL
  TYPE(const_mvec_ptr)  V;      // B-orthonormal basis
  TYPE(const_mvec_ptr)  BV;     // B*V
  TYPE(const_sdMat_ptr) sigma;  // matrix of shifts (diagonal for complex JDQR, block-diagonal for real JDQR with complex-conjugate ev, possibly upper triangular for block JDQR)
  TYPE(sdMat_ptr)       vy;     // temporary storage for V' times some block vector Y
  TYPE(mvec_ptr)        Work;   // temporary storage for block vectors, only needed for B != NULL
} TYPE(jadaOp_data);


// actually applies the jada-operator: X -> Y = (I-BV*V')(A X - B*X*sigma)
void SUBR(jadaOp_apply)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
{
  ENTER_FCN(__FUNCTION__);

  // convert op to jadaOp_data
  const TYPE(jadaOp_data) *jadaOp = (const TYPE(jadaOp_data)*) op;

#if PHIST_OUTLEV>=PHIST_DEBUG

int nrX,ncX,nrY,ncY,nrV,ncV,nrSig,ncSig;
PHIST_CHK_IERR(SUBR(mvec_my_length)(X,&nrX,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&ncX,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_my_length)(Y,&nrY,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Y,&ncY,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_my_length)(jadaOp->V,&nrV,ierr),*ierr);
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(jadaOp->V,&ncV,ierr),*ierr);

PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(jadaOp->sigma,&nrSig,ierr),*ierr);
PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(jadaOp->sigma,&ncSig,ierr),*ierr);

PHIST_DEB("X is %dx%d", nrX, ncX);
PHIST_DEB("Y is %dx%d", nrY, ncY);
PHIST_DEB("V is %dx%d", nrV, ncV);
PHIST_DEB("sigma is %dx%d", nrSig, ncSig);

#endif

  // Y = A*X - Y

  // first calculate Y = A*X - B*X*sigma
  // apply B if we have one
  if( jadaOp->B_op != NULL )
  {
    PHIST_CHK_IERR(jadaOp->B_op->apply(ONE,jadaOp->B_op->A,X,ZERO,jadaOp->Work,ierr),*ierr);      // Work = B*X
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(ONE,jadaOp->Work,jadaOp->sigma,ZERO,Y,ierr),*ierr);     // Y    = Work*sigma
  }
  else
  {
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(ONE,X,jadaOp->sigma,ZERO,Y,ierr),*ierr);                // Y    = X*sigma
  }

  PHIST_CHK_IERR(jadaOp->A_op->apply(ONE,jadaOp->A_op->A,X,-ONE,Y,ierr),*ierr);                   // Y = A*X - Y


  // then calculate Y = (Y - BV*V'*Y)
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(ONE,jadaOp->V,Y,ZERO,jadaOp->vy,ierr),*ierr);             // vy = V'*Y
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-ONE,jadaOp->BV,jadaOp->vy,ONE,Y,ierr),*ierr);            // Y  = Y - BV*vy
  return;
}


// allocate and initialize the jadaOp struct
void SUBR(jadaOp_create)(TYPE(const_op_ptr)    A_op,  TYPE(const_op_ptr)   B_op,
                         TYPE(const_mvec_ptr)  V,     TYPE(const_mvec_ptr) BV,
                         TYPE(const_sdMat_ptr) sigma, TYPE(mvec_ptr)       Work,
                         TYPE(op_ptr)*         jdOp,  int*                 ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;

  // allocate jadaOp struct
  *jdOp = (TYPE(op_ptr))malloc(sizeof(TYPE(op)));
  TYPE(jadaOp_data) *myOp = (TYPE(jadaOp_data)*)malloc(sizeof(TYPE(jadaOp_data)));

  // setup jadaOp members
  myOp->A_op  = A_op;

  myOp->B_op  = B_op;
  myOp->V     = V;
  myOp->sigma = sigma;

  // if B_op == NULL we can use V here
  if( B_op != NULL )
  {
    myOp->BV    = BV;
    myOp->Work  = Work;
  }
  else
  {
    myOp->BV    = V;
    myOp->Work  = NULL;
  }

  // allocate data for temporary sdMat vy
  int nvec_V, nvec_X;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec_V,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(sigma,&nvec_X,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&(myOp->vy),nvec_V,nvec_X,NULL,ierr),*ierr);

  // setup op_ptr function pointers
  (*jdOp)->A     = (const void*)myOp;
// note: this will cause a segfault later on if we don't explicitly cast to the correct 
// function return type (void*). This is needed here because we're in C, not C++.
//(*jdOp)->apply = (void(*)(_ST_,const void*,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*))
(*jdOp)->apply = (&SUBR(jadaOp_apply));

  return;
}


// deallocate jadaOp struct
void SUBR(jaraOp_delete)(TYPE(op_ptr)* jdOp, int *ierr)
{
  ENTER_FCN(__FUNCTION__);

  // get jadaOp
  TYPE(jadaOp_data) *jadaOp = (TYPE(jadaOp_data)*) jdOp;

  // delete op
  free(jdOp);
  jdOp = NULL;

  // delete temporary sdMat vy
  PHIST_CHK_IERR(SUBR(sdMat_delete)(jadaOp->vy,ierr),*ierr);

  // delete jadaOp
  free(jadaOp);

  return;
}

