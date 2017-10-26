/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
class TYPE(jadaOp_data)
{
  public:
  
  TYPE(jadaOp_data)(){}

  TYPE(const_linearOp_ptr)    AB_op;   // operator of the general matrix A
  TYPE(const_linearOp_ptr)    B_op;   // operator of the hpd. matrix B, assumed I when NULL
  TYPE(const_linearOp_ptr)    leftPrecon_op;   // left preconditioning operator
  TYPE(const_mvec_ptr)  V;      // B-orthonormal basis
  TYPE(const_mvec_ptr)  BV;     // B*V
  int num_shifts;               // number of shifts given to constructor
  const _ST_*           sigma;  // array of NEGATIVE shifts, assumed to have correct size; TODO: what about 'complex' shifts for real JDQR?
  TYPE(mvec_ptr)        X_proj; // temporary storage for (I-VV'B)X, only used for B!= NULL
  MvecOwner<_ST_> _X_proj, _V_prec; // these objects make sure that some temporary storage is freed when the object is deleted
};


//! simply apply original operator shifted, no pre- or postprojection
void SUBR(jadaOp_apply_project_none)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  int nvecX;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvecX,iflag),*iflag);
  if (nvecX>jadaOp->num_shifts)
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }

  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
  }
  else
  {

    // y_i <- alpha*(A+sigma_i I)*x_i + beta * y_i
    if (jadaOp->leftPrecon_op==NULL)
    {
      PHIST_CHK_IERR(SUBR(linearOp_apply_shifted)(alpha, jadaOp->AB_op, jadaOp->sigma, X, beta, Y, iflag),*iflag);
    }
    else
    {
      PHIST_CHK_IERR(*iflag = (alpha==st::one() && beta==st::zero())?0: PHIST_NOT_IMPLEMENTED,*iflag);
      TYPE(mvec_ptr) opX=NULL;
      PHIST_CHK_IERR(SUBR(mvec_clone_shape)(&opX,X,iflag),*iflag);
      MvecOwner<_ST_> _opX(opX);
      PHIST_CHK_IERR(SUBR(linearOp_apply_shifted)(alpha, jadaOp->AB_op, jadaOp->sigma, X, 
                                                  beta, opX, iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_apply)(alpha, jadaOp->leftPrecon_op, opX, beta, Y, iflag),*iflag);
    }
  }
}

// applies the jada-operator as in apply_project_pre_post but skips the final projection
//
// for B==NULL:  Y <- alpha*(AX_ +  X_*sigma) + beta*Y
//               with X_=(I-VV')X
//
// for B!=NULL:  Y <- alpha*(AX_ + BX_*sigma) + beta*Y
//               with X_=(I-V(VB)')X
//
void SUBR(jadaOp_apply_project_pre)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  int nvecX;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvecX,iflag),*iflag);
  if (nvecX>jadaOp->num_shifts)
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }

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

    // apply shifted A and preconditioner
    PHIST_CHK_IERR(SUBR(jadaOp_apply_project_none)(alpha, op, X_proj,beta, Y, iflag), *iflag);
    if (X_proj && X_proj != jadaOp->X_proj)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(X_proj,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp,iflag),*iflag);
  }
}


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

  int nvecX;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvecX,iflag),*iflag);
  if (nvecX>jadaOp->num_shifts)
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }

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
    PHIST_CHK_IERR(SUBR(linearOp_apply_shifted)(alpha, jadaOp->AB_op, jadaOp->sigma, X, beta, Y, iflag),*iflag);
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

  int nvecX;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvecX,iflag),*iflag);
  if (nvecX>jadaOp->num_shifts)
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
  }

  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale       ) (Y, beta, iflag), *iflag);
    return;
  }

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

  MvecOwner<_ST_> _X_proj( (X_proj==jadaOp->X_proj)?NULL: X_proj);
  SdMatOwner<_ST_> _tmp(tmp);

  // for now we have BV, but in the long run I think we should get rid of that extra vector block (cf. #

  // using a fused kernel here would save some data traffic, but the corresponding kernel doesn't exist right now:
  // X_proj = X - V*(BV'X)

  // tmp <- (BV)'*X and X_proj = X
{
  PHIST_ENTER_FCN("phist_jadaOp_mvecT_times_mvec_and_copy_x");
  PHIST_CHK_IERR( SUBR( mvecT_times_mvec ) (st::one(),  jadaOp->BV,  X,   st::zero(), tmp, iflag), *iflag);
  PHIST_CHK_IERR( SUBR( mvec_add_mvec ) (st::one(), X, st::zero(), X_proj, iflag), *iflag);
}
    // X_proj <- X - V*tm
{
  PHIST_ENTER_FCN("phist_jadaOp_mvec_times_sdMat");
  PHIST_CHK_IERR( SUBR( mvec_times_sdMat ) (-st::one(), jadaOp->V, tmp, st::one(),  X_proj,   iflag), *iflag);
}

  // apply shifted A and post-project
  PHIST_CHK_IERR(SUBR(jadaOp_apply_project_post)(alpha, op, X_proj,beta, Y, iflag), *iflag);

  // note: the preconditioner may carry its own post-projection, too
  if (jadaOp->leftPrecon_op!=NULL)
  {
      PHIST_CHK_IERR(*iflag = (alpha==st::one() && beta==st::zero())?0: PHIST_NOT_IMPLEMENTED,*iflag);
      PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),Y,st::zero(),X_proj,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_apply)(st::one(), jadaOp->leftPrecon_op, X_proj, st::zero(), Y, iflag),*iflag);
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
  TYPE(jadaOp_data) *myOp = new TYPE(jadaOp_data);
  
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
  myOp->leftPrecon_op = NULL;
  myOp->V      = V;
  myOp->BV     = (BV != NULL ? BV     : V);
  myOp->num_shifts = nvec;
  myOp->sigma  = sigma;
  // allocate necessary temporary arrays
  int nvecp=0;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_comm(AB_op->domain_map, &comm, iflag), *iflag);
  if (V!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V, &nvecp, iflag), *iflag);
  }
  // this vector is created for a fixed number of vectors nvec to which the operator
  // is applied. if another blocksize occurs, temporary storage is used in the apply
  // function, which may mean some overhead!
  myOp->X_proj = NULL;

  // setup op_ptr function pointers. For standard EVP, just post-project. For generalized EVP,
  // pre- and postproject. X_proj is only needed if B!=NULL, otherwise we don't do any pre-
  // projection.
  jdOp->A     = (const void*)myOp;

  if (nvecp==0)
  {
    jdOp->apply = SUBR(jadaOp_apply_project_none);
  }
  else
  {
    if (B_op!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_create)(&myOp->X_proj,AB_op->domain_map,nvec,iflag),*iflag);
      // make sure X_proj gets deleted automatically when myOp is deleted
      myOp->_X_proj.set(myOp->X_proj);
      // if the user passes in a B operator and projection vectors V, he also has to provide BV
      PHIST_CHK_IERR(*iflag= (V!=NULL && BV==V)? PHIST_INVALID_INPUT: 0, *iflag);
      jdOp->apply = SUBR(jadaOp_apply_project_pre_post);
    }
    else
    {
      jdOp->apply = SUBR(jadaOp_apply_project_post);
    }
  }
  jdOp->applyT= NULL; // not needed, I think, but it's trivial to implement
  jdOp->apply_shifted=NULL;// does not make sense, it would mean calling apply_shifted in a 
                           // nested way.
  jdOp->destroy = SUBR(jadaOp_delete);

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

  if (jdOp->A == NULL) return;

  // get jadaOp_data
  TYPE(jadaOp_data) *jadaOp = (TYPE(jadaOp_data)*) jdOp->A;

  // delete jadaOp
  delete jadaOp;
}

// create a preconditioner for the inner solve in Jacobi-Davidson.
//
// Given a linear operator that is a preconditioner for A-sigma_j*B, this function will
// wrap it up to use apply_shifted and pre-/post-apply adequate (skew-) projections
// when apply() is called. We need this because our implementations     
// of blockedGMRES and MINRES are not aware of the shifts so they can only call apply in the precon-    
// ditioning operator. Obviously not all preconditioners are able to handle varying shifts without      
// recomputing, this is not taken into account by this function:in that case the input P_op must be     
// updated beforehand, or the existing preconditioner for e.g. A or A-tau*B is applied.
//
// If V is given, the preconditioner application will include a skew-projection
//
// Y <- (I - P_op\V (BV'P_op\V)^{-1} (BV)') P_op\X
//
// If BV==NULL, BV=V is assumed.
extern "C" void SUBR(jadaPrec_create)(TYPE(const_linearOp_ptr) P_op,
        TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) BV,
        const _ST_ sigma[], int nvec,
        TYPE(linearOp_ptr) jdPrec, int projType, int* iflag)

{
#include "phist_std_typedefs.hpp"
  if (V==NULL||projType==0)
  {
    // simply apply the preconditioner "as is" or apply given projection vectors after applying P
    PHIST_CHK_IERR(SUBR(jadaOp_create)(P_op,NULL,V,BV,sigma, nvec, jdPrec,iflag),*iflag);
    // use the version without the projections:
    jdPrec->apply=SUBR(jadaOp_apply_project_none);
  }
  else
  {
    // we always require BV!=NULL, but it may point to the same memory or object in memory
    if (BV==NULL)
    {
      PHIST_SOUT(PHIST_ERROR,"jadaPrec_create requires BV if V is given, but it may be the same object or a view \n"
                             " of the same memory as V\n");
      *iflag=PHIST_INVALID_INPUT;
      return;
    }
    // construct a skew-projected operator, first we need to construct P\V*(BV'P\V)^{-1}
    int nproj;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nproj,iflag),*iflag);
    if (nproj==0) PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
    
    TYPE(mvec_ptr) PV=NULL;
    PHIST_CHK_IERR(SUBR(mvec_create)(&PV,P_op->domain_map,nproj,iflag),*iflag);
    PHIST_CHK_IERR(P_op->apply(st::one(),P_op->A,V,st::zero(),PV,iflag),*iflag);
    TYPE(sdMat_ptr) BVtPV=NULL, M=NULL;
    phist_const_comm_ptr comm=NULL;

    PHIST_CHK_IERR(phist_map_get_comm(P_op->domain_map,&comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&BVtPV, nproj,nproj,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&M, nproj,nproj,comm,iflag),*iflag);
    SdMatOwner<_ST_> _BVtPV(BVtPV), _M(M);
    
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),BV,PV,st::zero(),BVtPV,iflag),*iflag);

    // compute the *transposed* pseudo-inverse of V'K\V in place
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(BVtPV,iflag),*iflag);    
    int rank;
    _MT_ rankTol=mt::eps();
    PHIST_CHK_IERR(SUBR(sdMat_pseudo_inverse)(BVtPV,&rank,rankTol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(BVtPV,iflag),*iflag);
    // explicitly transpose the result
    PHIST_CHK_IERR(SUBR(sdMatT_add_sdMat)(st::one(),BVtPV,st::zero(),M,iflag),*iflag);

    // in-place PV*(V'P\V)^+
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(PV,M,iflag),*iflag);

    // use the version with only post-projection by giving B_op==NULL:
    PHIST_CHK_IERR(SUBR(jadaOp_create)(P_op,NULL,BV,PV,sigma, nvec, jdPrec,iflag),*iflag);
    jdPrec->apply=SUBR(jadaOp_apply_project_post);
    // add the PV block vector to be deleted automatically when the operator is deleted.
    TYPE(jadaOp_data)* jdDat=(TYPE(jadaOp_data)*)jdPrec->A;
    jdDat->_V_prec.set(PV);
  }
}

//! add a left preconditioner created by jadaPrec_create to a jadaOp.

//! The effect of the apply function will afterwards by Y <- alpha*(jadaPrec*jadaOp*X) + beta*Y,
//! the projections used are determined by the AB_op and jadaPrec operators. If jadaPrec==NULL, 
//! the operator is reset to it's original effect.
extern "C" void SUBR(jadaOp_set_leftPrecond)(TYPE(linearOp_ptr) jdOp, TYPE(const_linearOp_ptr) jadaPrec, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(TYPE(jadaOp_data), jdDat, jdOp->A, *iflag);
  jdDat->leftPrecon_op = jadaPrec;
  if (jadaPrec==NULL)
  {
    // recover original behavior: post-project for B=I, pre/post for B!=I
    if (jdDat->B_op==NULL)
    {
      jdOp->apply=SUBR(jadaOp_apply_project_post);
    }
    else
    {
      jdOp->apply=SUBR(jadaOp_apply_project_pre_post);
    }
  }
  else
  { 
    PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jdPrec, jadaPrec->A, *iflag);
    // We pre-project with (I-(BV)V') and
    // postproject with (I-Z(Z'BQ)^{-1}BQ^T), Z=P\Q. The post-projection
    // is handled by the preconditioner wrapper (see jadaPrec_create)
    jdOp->apply=SUBR(jadaOp_apply_project_pre_post);
  }
}
