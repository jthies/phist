
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
//! Remark: for the generalized eigenvalue problem with hpd. B
//!         we need always both projections, e.g. Y = (I-BVV')(A*X_-B*X_*sigma)
//!         with X_ = (I-VV'B)X
//!         (melven)
//! Remark: if we directly store -sigma, we can avoid an additional vector operation
//!         (as the mvec_vadd_mvec doesn't provei an additional scaling argument)


// private struct to keep all the pointers we need in order to apply the operator.
typedef struct TYPE(jadaOp_data)
{
  TYPE(const_linearOp_ptr)    A_op;   // operator of the general matrix A
  TYPE(const_linearOp_ptr)    B_op;   // operator of the hpd. matrix B, assumed I when NULL
  TYPE(const_mvec_ptr)  V;      // B-orthonormal basis
  TYPE(const_mvec_ptr)  BV;     // B*V
  const _ST_*           sigma;  // array of NEGATIVE shifts, assumed to have correct size; TODO: what about 'complex' shifts for real JDQR?
  TYPE(mvec_ptr)        X_proj; // temporary storage for (I-VV'B)X, only used for B!= NULL
} TYPE(jadaOp_data);


// actually applies the jada-operator:
// for B==NULL:  Y <- alpha* (I-VV')(AX+BX*sigma) + beta*Y
// for B!=NULL:  Y <- alpha* (I-BVV')*(A(I-VV'B)X + BX*sigma) + beta*Y
void SUBR(jadaOp_apply)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  PHIST_CHK_IERR(*iflag = (jadaOp->B_op != NULL) ? PHIST_NOT_IMPLEMENTED : 0, *iflag);

  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
  }
  else
  {
    int nvec, nvecp;
    PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (X,          &nvec,  iflag), *iflag);
    PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (jadaOp->V,  &nvecp, iflag), *iflag);
    const_comm_ptr_t comm;
    PHIST_CHK_IERR( SUBR( mvec_get_comm ) (X, &comm, iflag), *iflag);
    TYPE(sdMat_ptr) tmp;
    PHIST_CHK_IERR( SUBR( sdMat_create ) (&tmp, nvecp, nvec, comm, iflag), *iflag);

    // y_i <- alpha*(A+sigma_i I)*x_i + beta * y_i
{
PHIST_ENTER_FCN("phist_jadaOp_shifted_A_times_mvec");
    PHIST_CHK_IERR(jadaOp->A_op->apply_shifted(alpha, jadaOp->A_op->A, jadaOp->sigma, X, beta, Y, iflag),*iflag);
}
    // tmp <- V'*Y
{
PHIST_ENTER_FCN("phist_jadaOp_mvecT_times_mvec");
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(),  jadaOp->V,  Y,   st::zero(), tmp, iflag), *iflag);
}
    // Y <- Y - V*tmp
{
PHIST_ENTER_FCN("phist_jadaOp_mvec_times_sdMat");
    PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), jadaOp->BV, tmp, st::one(),  Y,   iflag), *iflag);
}
    PHIST_CHK_IERR( SUBR( sdMat_delete ) (tmp, iflag), *iflag);
  }

/*
  TYPE(const_mvec_ptr) BX;
  if( jadaOp->B_op == NULL )
  {
    // calculate AX and set BX = 
    PHIST_CHK_IERR( jadaOp->A_op->apply(st::one(), jadaOp->A_op->A, X, st::zero(), jadaOp->AX, iflag), *iflag);               // AX     <- A*X
    BX = X;
  }
  else // B_op != NULL
  {
    // calculate X_proj, AX_proj and BX_proj
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(), jadaOp->BV, X, st::zero(), jadaOp->VY, iflag), *iflag);              // VY     <- (BV)'X
    PHIST_CHK_IERR( SUBR( mvec_add_mvec    ) (st::one(), X, st::zero(), jadaOp->X_proj, iflag), *iflag);                      // X_proj <- X
    PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), jadaOp->V, jadaOp->VY, st::one(), jadaOp->X_proj, iflag), *iflag);  // X_proj <- X_proj - V*VY
    PHIST_CHK_IERR( jadaOp->A_op->apply(st::one(), jadaOp->A_op->A, jadaOp->X_proj, st::zero(), jadaOp->AX, iflag), *iflag);  // AX     <- A*X_proj
    PHIST_CHK_IERR( jadaOp->A_op->apply(st::one(), jadaOp->B_op->A, jadaOp->X_proj, st::zero(), jadaOp->BX, iflag), *iflag);  // BX     <- B*X_proj
    BX = jadaOp->BX;
  }


  if( alpha != st::zero() )
  {
    // Y <- alpha* (I-BVV')(AX-BX*sigma) + beta*Y, assumes (I-BVV')Y = Y
    PHIST_CHK_IERR( SUBR( mvec_add_mvec    ) (st::one(), jadaOp->AX, beta/alpha, Y, iflag), *iflag);                    // Y      <- AX + beta/alpha*Y
    PHIST_CHK_IERR( SUBR( mvec_vadd_mvec   ) (jadaOp->sigma, BX, st::one(), Y, iflag), *iflag);                         // Y      <- Y + BX*sigma
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(), jadaOp->V, Y, st::zero(), jadaOp->VY, iflag), *iflag);               // VY     <- V'*Y
    PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-alpha, jadaOp->BV, jadaOp->VY, alpha, Y, iflag), *iflag);          // Y      <- alpha*(Y - BV*VY)
  }
  else
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
  }
*/
}


// allocate and initialize the jadaOp struct
void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    A_op,    TYPE(const_linearOp_ptr)    B_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], int                   nvec,
                         TYPE(linearOp_ptr)          jdOp,    int*                  iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  int i;
  // allocate jadaOp struct
  TYPE(jadaOp_data) *myOp = new(TYPE(jadaOp_data));
  
  if (A_op->apply_shifted==NULL)
  {
    PHIST_SOUT(PHIST_ERROR, "operator passed to %s does not support apply_shifted\n"
                            "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
  *iflag=-1;
  return;
  }

  // setup jadaOp members
  myOp->A_op   = A_op;
  myOp->B_op   = B_op;
  myOp->V      = V;
  myOp->BV     = (B_op != NULL ? BV     : V);
  myOp->sigma  = sigma;
  // allocate necessary temporary arrays
  int nvecp;
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(phist_map_get_comm(A_op->domain_map, &comm, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V, &nvecp, iflag), *iflag);
  if( B_op != NULL )
  {
    PHIST_CHK_IERR(SUBR(mvec_create)(&myOp->X_proj, B_op->domain_map, nvec, iflag), *iflag);
  }
  else
  {
    myOp->X_proj = NULL;
  }

  // setup op_ptr function pointers
  jdOp->A     = (const void*)myOp;
  jdOp->apply = (&SUBR(jadaOp_apply));
  jdOp->applyT= NULL; // not needed, I think, but it's trivial to implement
  jdOp->apply_shifted=NULL;// does not make sense, it would mean calling apply_shifted in a 
                           // nested way.

  // print some useful data
  PHIST_SOUT(PHIST_DEBUG, "Created jadaOp with %d projection vectors and %d shifts\n",   nvecp,nvec);
}


// deallocate jadaOp struct
void SUBR(jadaOp_delete)(TYPE(linearOp_ptr) jdOp, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  if( jdOp == NULL )
    return;

  // get jadaOp
  TYPE(jadaOp_data) *jadaOp = (TYPE(jadaOp_data)*) jdOp->A;
  if( jdOp->A == NULL )
    return;

  // delete temporary arrays
  PHIST_CHK_IERR(SUBR(mvec_delete)(jadaOp->X_proj, iflag), *iflag);

  // delete jadaOp
  delete jadaOp;
}

