// this file implements various forms of (skew-)projected       
// linear operators that can be used in Jacobi-Davidson type    
// algorithms. The general form is:                             
//                                                              
// post-project:   y <- (I - W V')*Op          * x  (1)         
// pre-/post:      y <- (I - W V')*Op*(I-V W') * x  (2)         
//                                                              
// In Jacobi-Davidson, we want to perform an iterative solve    
// (e.g. GMRES) for the correction equation. The                
// corresponding operator we put in the GMRES solver is         
//                                                              
// (I-QQ')(A-sigma_j I)               for standard EVP, and     
// (I-(BQ)Q')(A-sigma_j B)(I-Q(BQ)')  for generalized EVP with  
//                              hpd B. Here Q'BQ=I.             
//                                                              
// Q = [Q_locked v] is the basis of locked eigenvectors,        
//  extended by the current search direction(s) v (also called  
// Qtil in the JaDa implementations).                           
//                                                              
// We want the GMRES method to run for several shifts           
// sigma_j simultaneously, which we achieve by using the        
// operator's apply_shifted function.                           
//                                                              
// skew-projected preconditioning                               
// ==============================                               
//                                                              
/// If the jada option "preconSkewProject" is not 0, the correc-
// tion equation is preconditioned using (here K is a precondi- 
// tioner for (A-sigma*B):                                      
//                                                              
// x <- (I - (K\V)*((BV)'K\V)^{-1} (BV)') K\y   (3)             
//                                                              
// To save memory and projection operations we use V=v instead  
// of V=Q (which would project all converged eigenvectors out). 
// This can be implemented by (1) using W = (K\v)((Bv)'K\v)^{-1}
// and V=Bv. If K is symmetric and we want to preserve this     
// property, pre-/post (2) could be used, too.                  
//                                                              
// derivation:                                                  
//                                                              
// correction equation with B=B' hpd:                           
//                                                              
// (A - sigma*B)t = Q*S -r      with S s.t.    Q'Bt=0           
//                                                              
// approximate LHS by a preconditioner K:                       
//                                                              
// t_      = K \ (Q*S - r_)            s.t. t_'Bv=0             
//                                                              
// Q'Bt_   = (Q'BK\Q)*S - Q'BK\r_      = 0                      
// =>    S = (Q'B*K\Q)^{-1} (Q'BK\r_)                           
//                                                              
// Let y:=-r_ and x:=t_ here.                                   
// this results in the preconditioning operator (3)             
//                                                              



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
extern "C" void SUBR(jadaOp_create)(TYPE(const_linearOp_ptr)    AB_op,
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
extern "C" void SUBR(jadaOp_delete)(TYPE(linearOp_ptr) jdOp, int *iflag)
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

// create a preconditioner for the inner solve in Jacobi-Davidson.
//
// Given a linear operator that is a preconditioner for A, this function will simply
// wrap it up to use apply_shifted when apply() is called. We need this because our implementations
// of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-
// ditioning operator. Obviously not all preconditioners are able to handle varying shifts without
// recomputing, this is not taken into account by this function:in that case the input P_op must be
// updated beforehand.
//
// If V is given, the preconditioner application will include a skew-projection
// Y <- (I - P_op\V (V' P_op\ V)^{-1} V' ) P_op\X      or (if BV!=V and BV!=NULL):
// Y <- (I -BP_op\V (V'BP_op\BV)^{-1} V'B) P_op\X
//
extern "C" void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op, 
        TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) BV,
        const _ST_ sigma[], int nvec,
        TYPE(linearOp_ptr) jdPrec, int* iflag)
{
#include "phist_std_typedefs.hpp"
  if (V==NULL)
  {
    // simply apply the preconditioner "as is"
    PHIST_CHK_IERR(SUBR(jadaOp_create)(P_op,NULL,NULL,NULL,sigma, nvec, jdPrec,iflag),*iflag);
    // use the version without the projections:
    jdPrec->apply=SUBR(jadaOp_apply_project_none);
  }
  else
  {
    // construct a skew-projected operator, first we need to construct P\V*(V'P\V)^{-1}
    if (BV!=V && BV!=NULL)
    {
      PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
    }
    else
    {
      int nproj;
      PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nproj,iflag),*iflag);
      if (nproj==0) PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
      
      TYPE(mvec_ptr) PV=NULL;
      PHIST_CHK_IERR(SUBR(mvec_create)(&PV,P_op->domain_map,nproj,iflag),*iflag);
      PHIST_CHK_IERR(P_op->apply(st::one(),P_op->A,V,st::zero(),PV,iflag),*iflag);
      TYPE(sdMat_ptr) VtPV=NULL;
      phist_const_comm_ptr comm=NULL;
      PHIST_CHK_IERR(phist_map_get_comm(P_op->domain_map,&comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&VtPV, nproj,nproj,comm,iflag),*iflag);
      
      // compute the pseudo-inverse of V'K\V in place
      int rank;
      PHIST_CHK_IERR(SUBR(sdMat_pseudo_inverse)(VtPV,&rank,iflag),*iflag);
      
      // in-place PV*(V'P\V)^+
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(PV,VtPV,iflag),*iflag);
      
      // delete temporary sdMat
      PHIST_CHK_IERR(SUBR(sdMat_delete)(VtPV,iflag),*iflag);
      
      //TODO: when deleting this operator, PV must be deleted!
      //      We should introduce std::shared_ptr objects (available in C++11, it seems!)

      // use the version with only post-projection:
      PHIST_CHK_IERR(SUBR(jadaOp_create)(P_op,NULL,V,PV,sigma, nvec, jdPrec,iflag),*iflag);
      jdPrec->apply=SUBR(jadaOp_apply_project_post);
    }
  }
}
