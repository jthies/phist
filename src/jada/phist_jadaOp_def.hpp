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

  TYPE(linearOp_ptr)    AB_op;   // operator of the general matrix A
  TYPE(const_linearOp_ptr)    B_op;   // operator of the hpd. matrix B, assumed I when NULL
  TYPE(const_linearOp_ptr)    leftPrecon_op;   // left preconditioning operator
  TYPE(linearOp_ptr)    preProj_op;   // pre-projection operator
  TYPE(linearOp_ptr)    postProj_op;    // post-projection operator
  TYPE(linearOp_ptr)    Prec_op;   // projection operator
  TYPE(linearOp_ptr)    Skew_op;   // projection operator
  TYPE(linearOp_ptr)    k_op; //jadaOp with k_wrapper
  TYPE(const_mvec_ptr)  V;      // B-orthonormal basis
  TYPE(const_mvec_ptr)  BV;     // B*V
  int num_shifts;               // number of shifts given to constructor
  int projType;
  const _ST_*           sigma;  // array of NEGATIVE shifts, assumed to have correct size; TODO: what about 'complex' shifts for real JDQR?
  const _ST_*           sigma_prec; //array of shifts for the preconditioner Prec_op
};

// private struct to keep all the pointers we need in order to apply the operator.
class TYPE(projOp_data)
{
  public:

  TYPE(projOp_data)(){}

  TYPE(const_mvec_ptr)  V;
  TYPE(const_mvec_ptr)  W;
  MvecOwner<_ST_> _V, _W; // these objects make sure that some temporary storage is freed when the object is deleted
};

//! apply for jadaOp
void SUBR(jadaOp_apply)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X,
    _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const TYPE(jadaOp_data), jadaOp, op, *iflag);

  TYPE(const_linearOp_ptr) k_Op = jadaOp->k_op;

  if(k_Op == NULL)
  {
    PHIST_CHK_IERR(*iflag=PHIST_INVALID_INPUT,*iflag);
    return;
  }
  if( alpha == st::zero() )
  {
    PHIST_CHK_IERR( SUBR( mvec_scale ) (Y, beta, iflag), *iflag);
  }
  else
  {
    k_Op->apply(alpha,k_Op->A,X,beta,Y,iflag);
  }
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

  // delete inner Operators
  if (jadaOp->k_op!=NULL) PHIST_CHK_IERR(jadaOp->k_op->destroy(jadaOp->k_op, iflag), *iflag);
  delete jadaOp->k_op;

  if (jadaOp->preProj_op!=NULL) PHIST_CHK_IERR(jadaOp->preProj_op->destroy(jadaOp->preProj_op, iflag), *iflag);
  delete jadaOp->preProj_op;
  
  if (jadaOp->postProj_op!=NULL) PHIST_CHK_IERR(jadaOp->postProj_op->destroy(jadaOp->postProj_op, iflag), *iflag);
  delete jadaOp->postProj_op;

  if (jadaOp->Skew_op!=NULL) PHIST_CHK_IERR(jadaOp->Skew_op->destroy(jadaOp->Skew_op, iflag), *iflag);
  delete jadaOp->Skew_op;
  
  delete jadaOp->AB_op;
  
  delete jadaOp->Prec_op;
  
  // delete jadaOp
  delete jadaOp;
}

// apply function for projection Operator
// Y <- alpha*(I-V*W')*X + beta*Y
void SUBR(projection_Op_apply)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    PHIST_CAST_PTR_FROM_VOID(const TYPE(projOp_data), proj_Op, op, *iflag);

    if( alpha == st::zero() ) {
        PHIST_CHK_IERR( SUBR( mvec_scale)(Y, beta, iflag), *iflag);
    }
    else {
        int nvec, nvecp = 0;
        PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (X, &nvec, iflag), *iflag);
        PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (proj_Op->W, &nvecp, iflag), *iflag);

        phist_const_map_ptr map=NULL;
        PHIST_CHK_IERR(SUBR(mvec_get_map)(X,&map,iflag),*iflag);
        phist_const_comm_ptr comm=NULL;
        PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);

        TYPE(sdMat_ptr) tmp;
        PHIST_CHK_IERR( SUBR( sdMat_create)(&tmp, nvecp, nvec, comm, iflag), *iflag);

        // tmp <- -W'*X
{
PHIST_ENTER_FCN("phist_jadaOp_mvecT_times_mvec");
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec)(-st::one(), proj_Op->W, X, st::zero(), tmp, iflag), *iflag);
}
        // Y <- alpha*X + beta*Y
        PHIST_CHK_IERR( SUBR(mvec_add_mvec)(alpha, X, beta, Y, iflag), *iflag);

        // Y <- Y + alpha*V*tmp
{
PHIST_ENTER_FCN("phist_jadaOp_mvec_times_sdMat");
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,proj_Op->V,tmp,st::one(),Y,iflag),*iflag);
}

        //clean up
        PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp,iflag),*iflag);
    }
}

// applyT function for projection Operator
// Y <- alpha*(I-W*V')*X + beta*Y
void SUBR(projection_Op_applyT)(_ST_ alpha, const void* op, TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    PHIST_CAST_PTR_FROM_VOID(const TYPE(projOp_data), proj_Op, op, *iflag);

    if( alpha == st::zero() ) {
        PHIST_CHK_IERR( SUBR( mvec_scale)(Y, beta, iflag), *iflag);
    }
    else {
        int nvec, nvecp = 0;
        PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (X, &nvec, iflag), *iflag);
        PHIST_CHK_IERR( SUBR( mvec_num_vectors ) (proj_Op->V, &nvecp, iflag), *iflag);

        phist_const_map_ptr map=NULL;
        PHIST_CHK_IERR(SUBR(mvec_get_map)(X,&map,iflag),*iflag);
        phist_const_comm_ptr comm=NULL;
        PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);

        TYPE(sdMat_ptr) tmp;
        PHIST_CHK_IERR( SUBR( sdMat_create)(&tmp, nvecp, nvec, comm, iflag), *iflag);

        // tmp <- -V'*X
{
PHIST_ENTER_FCN("phist_jadaOp_mvecT_times_mvec");
    PHIST_CHK_IERR( SUBR( mvecT_times_mvec)(-st::one(), proj_Op->V, X, st::zero(), tmp, iflag), *iflag);
}
        // Y <- alpha*X + beta*Y
        PHIST_CHK_IERR( SUBR(mvec_add_mvec)(alpha, X, beta, Y, iflag), *iflag);

        // Y <- Y + alpha*W*tmp
{
PHIST_ENTER_FCN("phist_jadaOp_mvec_times_sdMat");
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,proj_Op->W,tmp,st::one(),Y,iflag),*iflag);
}

        //clean up
        PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp,iflag),*iflag);
    }
}

// deallocate projOp struct
extern "C" void SUBR(projection_Op_delete)(TYPE(linearOp_ptr) proj_Op, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  if( proj_Op == NULL )
    return;

  if (proj_Op->A == NULL) return;

  // get jadaOp_data
  TYPE(projOp_data) *projOp = (TYPE(projOp_data)*) proj_Op->A;

  // delete jadaOp
  delete projOp;
}

// allocate and initialize the projOp struct
extern "C" void SUBR(projection_Op_create)(TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) W, TYPE(linearOp_ptr) proj_Op, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag = 0;

    // allocate projOp struct
    TYPE(projOp_data) *myOp = new TYPE(projOp_data);

    // setup jadaOp members
    myOp->V      = V;
    myOp->W     = W;

    proj_Op->A     = (const void*)myOp;
    proj_Op->apply = SUBR(projection_Op_apply);
    proj_Op->applyT = SUBR(projection_Op_applyT);
    proj_Op->apply_shifted = NULL;

    proj_Op->destroy = SUBR(projection_Op_delete);

    phist_const_map_ptr mapV=NULL;
    phist_const_map_ptr mapW=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&mapV,iflag),*iflag)
    PHIST_CHK_IERR(SUBR(mvec_get_map)(W,&mapW,iflag),*iflag)
    proj_Op->range_map = mapW;
    proj_Op->domain_map = mapV;
    proj_Op->use_transpose = 0;
    proj_Op->shifts = NULL;

}

extern "C" void SUBR(skew_projection_Op_create)(TYPE(const_linearOp_ptr) P_op,
        TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) BV,
        TYPE(linearOp_ptr) skew_Op, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag = 0;

    // we always require BV!=NULL, but it may point to the same memory or object in memory as V
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

    // compute the *transposed* pseudo-inverse of (BV)'P\V in place
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(BVtPV,iflag),*iflag);
    int rank;
    _MT_ rankTol=mt::eps();
    PHIST_CHK_IERR(SUBR(sdMat_pseudo_inverse)(BVtPV,&rank,rankTol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(BVtPV,iflag),*iflag);
    // explicitly transpose the result
    PHIST_CHK_IERR(SUBR(sdMatT_add_sdMat)(st::one(),BVtPV,st::zero(),M,iflag),*iflag);

    // in-place PV*((BV)'P\V)^+
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(PV,M,iflag),*iflag);

    PHIST_CHK_IERR(SUBR(projection_Op_create)(PV, BV, skew_Op, iflag),*iflag);

    // add the PV block vector to be deleted automatically when the operator is deleted.
    TYPE(projOp_data)* pjDat=(TYPE(projOp_data)*)skew_Op->A;
    pjDat->_W.set(PV);
}

// allocate and initialize the jadaOp struct for a variable combination of projections
extern "C" void SUBR(jadaOp_variable_create)(TYPE(const_linearOp_ptr)    AB_op,
                         TYPE(const_linearOp_ptr)     B_op, TYPE(const_linearOp_ptr)    Prec_op,
                         TYPE(const_mvec_ptr)  V,       TYPE(const_mvec_ptr)  BV,
                         const _ST_            sigma[], const _ST_ sigma_prec[], int    nvec,
                         TYPE(linearOp_ptr)          jdOp, phist_Eprojection method,
                         int onlyPrec, int*                  iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  // allocate jadaOp struct
  TYPE(jadaOp_data) *myOp = new TYPE(jadaOp_data); 

  if (AB_op->apply_shifted==NULL)
  {
    PHIST_SOUT(PHIST_ERROR, "operator passed to %s does not support apply_shifted\n"
                            "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
  *iflag=-1;
  return;
  }
  
  myOp->AB_op = new TYPE(linearOp);
  *(myOp->AB_op) = *AB_op;
  myOp->AB_op->shifts = (_ST_*)sigma;
  myOp->AB_op->use_transpose=0;
  
  if(Prec_op != NULL)
  {
    myOp->Prec_op = new TYPE(linearOp);
    *(myOp->Prec_op) = *Prec_op;
    myOp->Prec_op->shifts = (_ST_*)sigma_prec;
    myOp->Prec_op->use_transpose=0;
  }
  
  // setup jadaOp members
  myOp->B_op   = B_op;
  myOp->leftPrecon_op = NULL;
  myOp->V      = V;
  myOp->BV     = (BV != NULL ? BV     : V);
  myOp->num_shifts = nvec;
  myOp->sigma  = sigma;
  myOp->sigma_prec = sigma_prec;
  myOp->projType = onlyPrec;
  myOp->k_op = new TYPE(linearOp);
  PHIST_CHK_IERR(SUBR(linearOp_product_create)(myOp->k_op,iflag),*iflag);
  
  int nvecp=0;
  if (V!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V, &nvecp, iflag), *iflag);
  }  
  
  // case NONE
  if(method==phist_PROJ_NONE)
  {
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
  }
  
  // case PRE
  else if(method==phist_PROJ_PRE)
  {
    myOp->preProj_op = new TYPE(linearOp);
    PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->preProj_op,iflag),*iflag);
    
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->preProj_op,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
  }
  
  // case POST
  else if(method==phist_PROJ_POST)
  {
    myOp->postProj_op = new TYPE(linearOp);
    PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->postProj_op,iflag),*iflag);
    myOp->postProj_op->use_transpose = 1;
  
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->postProj_op,iflag),*iflag);
  }

  // case PRE_POST
  else if(method==phist_PROJ_PRE_POST)
  {
    myOp->preProj_op = new TYPE(linearOp);
    PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->preProj_op,iflag),*iflag);
    myOp->postProj_op = new TYPE(linearOp);
    PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->postProj_op,iflag),*iflag);
    myOp->postProj_op->use_transpose = 1;  
      
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->preProj_op,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->postProj_op,iflag),*iflag);
  }

  // case SKEW
  else if(method==phist_PROJ_SKEW)
  {
    if(onlyPrec == 1)
    {
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->Prec_op,iflag),*iflag);
    }
    else
    {
      myOp->Skew_op = new TYPE(linearOp);
      PHIST_CHK_IERR(SUBR(skew_projection_Op_create)(myOp->Prec_op,myOp->V,myOp->BV,myOp->Skew_op,iflag),*iflag);
    
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->Prec_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->Skew_op,iflag),*iflag);
    }
  }

  // case ALL
  else if(method==phist_PROJ_ALL)
  {
    if(onlyPrec == 1)
    {
      myOp->preProj_op = new TYPE(linearOp);
      PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->preProj_op,iflag),*iflag);
      myOp->postProj_op = new TYPE(linearOp);
      PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->postProj_op,iflag),*iflag);
      myOp->postProj_op->use_transpose = 1; 
    
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->preProj_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->postProj_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->Prec_op,iflag),*iflag);
    }
    else
    {
      myOp->preProj_op = new TYPE(linearOp);
      PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->preProj_op,iflag),*iflag);
      myOp->postProj_op = new TYPE(linearOp);
      PHIST_CHK_IERR(SUBR(projection_Op_create)(myOp->V,myOp->BV,myOp->postProj_op,iflag),*iflag);
      myOp->postProj_op->use_transpose = 1;
      myOp->Skew_op = new TYPE(linearOp);
      PHIST_CHK_IERR(SUBR(skew_projection_Op_create)(myOp->Prec_op,myOp->V,myOp->BV,myOp->Skew_op,iflag),*iflag);
        
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->preProj_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->AB_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->postProj_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->Prec_op,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(linearOp_product_extend)(myOp->k_op,myOp->Skew_op,iflag),*iflag);
    }
  }
  
  else
  {
    PHIST_SOUT(PHIST_ERROR,"there is no such method in %s implemented\n"
                           "choose between NONE 0, PRE 1, POST 2, PRE_POST 3, SKEW 4 and ALL 7\n"
                           "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
  *iflag=-1;
  return;
  } 

  jdOp->A = (const void*)myOp;
  jdOp->apply = SUBR(jadaOp_apply);
  jdOp->applyT = NULL;
  jdOp->apply_shifted = NULL;
  jdOp->destroy = SUBR(jadaOp_delete);

  jdOp->range_map = myOp->k_op->range_map;
  jdOp->domain_map = myOp->k_op->range_map;

  // print some useful data
  PHIST_SOUT(PHIST_DEBUG, "Created jadaOp with %d projection vectors and %d shifts\n%s\n",   nvecp,nvec,
              B_op==NULL? "                    B=I and postprojection\n":
                          "                    B-inner product and pre-/postprojection\n");  
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

  phist_Eprojection method;
  if (B_op!=NULL)
  {
    method = phist_PROJ_PRE_POST;
  }
  else
  {
    method = phist_PROJ_POST;
  }

  PHIST_CHK_IERR(SUBR(jadaOp_variable_create)(AB_op,B_op,NULL,V,BV,sigma,NULL,nvec,jdOp,method,0,iflag),*iflag);

}

extern "C" void SUBR(jadaPrec_delete)(TYPE(linearOp_ptr) jdPrec, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;

  if(jdPrec == NULL)
    return;

  if(jdPrec->A == NULL) return;

  //get jadaPrec_data
  TYPE(jadaOp_data) *jadaPrec = (TYPE(jadaOp_data)*) jdPrec->A;

  // delete inner Operators
  if (jadaPrec->k_op!=NULL) PHIST_CHK_IERR(jadaPrec->k_op->destroy(jadaPrec->k_op, iflag), *iflag);
  delete jadaPrec->k_op;

  if (jadaPrec->Skew_op!=NULL) PHIST_CHK_IERR(jadaPrec->Skew_op->destroy(jadaPrec->Skew_op, iflag), *iflag);
  delete jadaPrec->Skew_op;
  
  delete jadaPrec->AB_op;
  
  // delete jadaOp
  delete jadaPrec;
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
    // need methode NONE with P_op as AB_op
    phist_Eprojection method = phist_PROJ_NONE;
    PHIST_CHK_IERR(SUBR(jadaOp_variable_create)(P_op, NULL, NULL, V, BV, sigma, NULL, nvec, jdPrec, method, 1, iflag),*iflag);
    PHIST_CAST_PTR_FROM_VOID(TYPE(jadaOp_data), jdPr, jdPrec->A, *iflag);
    jdPr->projType=projType;
    jdPr->Prec_op = (TYPE(linearOp_ptr))P_op;
    jdPr->Skew_op = NULL;
    jdPrec->destroy = SUBR(jadaPrec_delete);
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

    phist_Eprojection method = phist_PROJ_NONE;
    PHIST_CHK_IERR(SUBR(jadaOp_variable_create)(P_op, NULL, NULL, V, BV, sigma, NULL, nvec, jdPrec, method, 0, iflag),*iflag);

    PHIST_CAST_PTR_FROM_VOID(TYPE(jadaOp_data), jdPr, jdPrec->A, *iflag);
    jdPr->projType=projType;
    jdPr->Prec_op = (TYPE(linearOp_ptr))P_op;
    jdPr->Skew_op = new TYPE(linearOp);
    PHIST_CHK_IERR(SUBR(skew_projection_Op_create)(P_op,V,BV,jdPr->Skew_op,iflag),*iflag);
    
    // jdPrec should be skew_op * p_op^shifted
    // from jadaOp_variable_create with method NONE we get p_op^shifted
    // so just need to append Skew_op
    PHIST_CHK_IERR(SUBR(linearOp_product_extend)(jdPr->k_op,jdPr->Skew_op,iflag),*iflag);
    jdPrec->destroy = SUBR(jadaPrec_delete);
  }
}

//! add a left preconditioner created by jadaPrec_create to a jadaOp.

//! The effect of the apply function will afterwards by Y <- alpha*(jadaPrec*jadaOp(with PRE_POST)*X) + beta*Y,
//! the projections used are determined by the AB_op and jadaPrec operators. If jadaPrec==NULL,
//! the operator is reset to it's original effect.
extern "C" void SUBR(jadaOp_set_leftPrecond)(TYPE(linearOp_ptr) jdOp, TYPE(const_linearOp_ptr) jadaPrec, int* iflag)
{
  if (jadaPrec==NULL)
  {
    return;
  }
  else
  {
    PHIST_CAST_PTR_FROM_VOID(TYPE(jadaOp_data), jdDat, jdOp->A, *iflag);
    PHIST_CAST_PTR_FROM_VOID(TYPE(jadaOp_data), jdPrec, jadaPrec->A, *iflag);  
    phist_Eprojection method= phist_PROJ_ALL;
    TYPE(linearOp) *jadaOp = new TYPE(linearOp);
    PHIST_CHK_IERR(SUBR(jadaOp_variable_create)(jdDat->AB_op, jdDat->B_op, jdPrec->Prec_op, jdDat->V, jdDat->BV, jdDat->sigma, jdPrec->sigma, jdDat->num_shifts, jadaOp, method, jdPrec->projType, iflag),*iflag);  
    PHIST_CAST_PTR_FROM_VOID(TYPE(jadaOp_data), jadaDat, jadaOp->A, *iflag);
    jadaDat->leftPrecon_op = jadaPrec;
    
    jdOp->destroy(jdOp,iflag);
    *jdOp = *jadaOp;
    delete jadaOp;
  }
}

