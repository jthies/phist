
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
  TYPE(const_linearOp_ptr)    AB_op;   // operator of the general matrix A
  TYPE(const_linearOp_ptr)    B_op;   // operator of the hpd. matrix B, assumed I when NULL
  TYPE(const_mvec_ptr)  V;      // B-orthonormal basis
  TYPE(const_mvec_ptr)  BV;     // B*V
  const _ST_*           sigma;  // array of NEGATIVE shifts, assumed to have correct size; TODO: what about 'complex' shifts for real JDQR?
  TYPE(mvec_ptr)        X_proj; // temporary storage for (I-VV'B)X, only used for B!= NULL
} TYPE(jadaOp_data);


// applies the jada-operator with only post-projection:
//
// for B==NULL:  Y <- alpha* (I-VV') * (AX +  X*sigma) + beta*Y
//
// for (B!=NULL: Y <- alpha* (I-BVV')* (AX + BX*sigma) + beta*Y
//
void SUBR(jadaOp_apply_project_post)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
  }
  else
  {
    int nvec, nvecp;
    PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (X,          &nvec,  iflag), *iflag);
    PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (jadaOp->V,  &nvecp, iflag), *iflag);
    phist_const_comm_ptr comm;
    PHIST_CHK_IERR( SUBR( mvec_get_comm ) (X, &comm, iflag), *iflag);
    TYPE(sdMat_ptr) tmp;
    PHIST_CHK_IERR( SUBR( sdMat_create ) (&tmp, nvecp, nvec, comm, iflag), *iflag);

    // y_i <- alpha*(A+sigma_i I)*x_i + beta * y_i
{
PHIST_ENTER_FCN("phist_jadaOp_shifted_A_times_mvec");
    PHIST_CHK_IERR(jadaOp->AB_op->apply_shifted(alpha, jadaOp->AB_op->A, jadaOp->sigma, X, beta, Y, iflag),*iflag);
}
    // tmp <- V'*Y
{
PHIST_ENTER_FCN("phist_jadaOp_mvecT_times_mvec");
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(),  jadaOp->V,  Y,   st::zero(), tmp, iflag), *iflag);
}
    // Y <- Y - BV*tmp
{
PHIST_ENTER_FCN("phist_jadaOp_mvec_times_sdMat");
    PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), jadaOp->BV, tmp, st::one(),  Y,   iflag), *iflag);
}
    PHIST_CHK_IERR( SUBR( sdMat_delete ) (tmp, iflag), *iflag);
  }
}

// applies the complete jada-operator:
//
// for B==NULL:  Y <- alpha* (I-VV') * (AX_ +  X_*sigma) + beta*Y
//               with X_=(I-VV')X
//
// for B!=NULL:  Y <- alpha* (I-BVV')* (AX_ + BX_*sigma) + beta*Y
//               with X_=(I-V(VB)')X
//
void SUBR(jadaOp_apply_project_pre_post)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
  }
  else
  {

    // pre-project X_proj <- (I-VV'B)X
    TYPE(mvec_ptr) X_proj = jadaOp->X_proj;
    int nvec, nvecp, nvec_tmp=0;
    PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (X,             &nvec,     iflag), *iflag);
    if (jadaOp->X_proj)
    {
      PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (jadaOp->X_proj, &nvec_tmp, iflag), *iflag);
    }
    PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (jadaOp->V,     &nvecp,    iflag), *iflag);

      phist_const_map_ptr map=NULL;
      PHIST_CHK_IERR(SUBR(mvec_get_map)(X,&map,iflag),*iflag);
    phist_const_comm_ptr comm=NULL;
    PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);
    
    if (nvec!=nvec_tmp || X_proj==NULL)
    {
      // can not use the temporary vector in the jadaOp
      PHIST_CHK_IERR(SUBR(mvec_create)(&X_proj,map,nvec,iflag),*iflag);      
    }
    
    TYPE(sdMat_ptr) tmp;
    PHIST_CHK_IERR( SUBR( sdMat_create ) (&tmp, nvecp, nvec, comm, iflag), *iflag);

    // for now we have BV, but in the long run I think we should get rid of that extra vector block (cf. #

    // using a fused kernel here would save some data traffic, but the corresponding kernel doesn't exist right now:
    // X_proj = X - V*(BV'X)

    // tmp <- (BV)'*X and X_proj = X
{
PHIST_ENTER_FCN("phist_jadaOp_mvecT_times_mvec_and_copy_x");
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(),  jadaOp->BV,  X,   st::zero(), tmp, iflag), *iflag);
    PHIST_CHK_IERR( SUBR( mvec_add_mvec ) (st::one(), X, st::zero(), X_proj, iflag), *iflag);
}
    // X_proj <- X - V*tmp
{
PHIST_ENTER_FCN("phist_jadaOp_mvec_times_sdMat");
    PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), jadaOp->V, tmp, st::one(),  X_proj,   iflag), 
    *iflag);
}

    // apply shifted A and post-project
    PHIST_CHK_IERR(SUBR(jadaOp_apply_project_post)(alpha, op, X_proj,beta, Y, iflag), *iflag);
    if (X_proj && X_proj != jadaOp->X_proj)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(X_proj,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp,iflag),*iflag);
  }
}


// allocate and initialize the jadaOp struct
void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op,
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
  
  if (AB_op->apply_shifted==NULL)
  {
    PHIST_SOUT(PHIST_ERROR, "operator passed to %s does not support apply_shifted\n"
                            "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
  *iflag=-1;
  return;
  }

  // setup jadaOp members
  myOp->AB_op   = AB_op;
  myOp->B_op   = B_op;
  myOp->V      = V;
  myOp->BV     = (BV != NULL ? BV     : V);
  myOp->sigma  = sigma;
  // allocate necessary temporary arrays
  int nvecp=0;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_comm(AB_op->domain_map, &comm, iflag), *iflag);
  if (V!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V, &nvecp, iflag), *iflag);
  }
  // this vector is created if needed in apply() (only if pre-projection is desired)
  myOp->X_proj = NULL;

  // setup op_ptr function pointers. For standard EVP, just post-project. For generalized EVP,
  // pre- and postproject.
  jdOp->A     = (const void*)myOp;
  if (B_op!=NULL)
  {
    // if the user passes in a B opeartor and projection vectors V, he also has to provide BV
    PHIST_CHK_IERR(*iflag= (V!=NULL && BV==V)? PHIST_INVALID_INPUT: 0, *iflag);
    jdOp->apply = (&SUBR(jadaOp_apply_project_pre_post));
  }
  else
  {
    jdOp->apply = (&SUBR(jadaOp_apply_project_post));
  }
  jdOp->applyT= NULL; // not needed, I think, but it's trivial to implement
  jdOp->apply_shifted=NULL;// does not make sense, it would mean calling apply_shifted in a 
                           // nested way.
  jdOp->destroy = (&SUBR(jadaOp_delete));

  jdOp->range_map=AB_op->range_map;
  jdOp->domain_map=AB_op->domain_map;
  
  // print some useful data
  PHIST_SOUT(PHIST_DEBUG, "Created jadaOp with %d projection vectors and %d shifts\n%s\n",   nvecp,nvec,
              B_op==NULL? "                    B=I and postprojection\n":
                          "                    B-inner product and pre-/postprojection\n");
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
  if (jadaOp->X_proj!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(jadaOp->X_proj, iflag), *iflag);
  }

  // delete jadaOp
  delete jadaOp;
}



void SUBR(jadaOp_apply_project_none)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
  }
  else
  {

    // y_i <- alpha*(A+sigma_i I)*x_i + beta * y_i
    PHIST_CHK_IERR(jadaOp->AB_op->apply_shifted(alpha, jadaOp->AB_op->A, jadaOp->sigma, X, beta, Y, iflag),*iflag);
  }
}

//! create a preconditioner for the inner solve in Jacobi-Davidson.
//!
//! Given a linear operator that is a preconditioner for A, this function will simply
//! wrap it up to use apply_shifted when apply() is called. We need this because our implementations
//! of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-
//! ditioning operator. Obviously not all preconditioners are able to handle varying shifts without
//! recomputing, this is not taken into account by this function:in that case the input P_op must be
//! updated beforehand.
void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op, const _ST_ sigma[], int nvec,
        TYPE(linearOp_ptr) jdPrec, int* iflag)
{
  PHIST_CHK_IERR(SUBR(jadaOp_create)(P_op,NULL,NULL,NULL,sigma, nvec, jdPrec,iflag),*iflag);
  // use the version without the projections:
  jdPrec->apply=SUBR(jadaOp_apply_project_none);
}
