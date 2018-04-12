/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
extern "C" {

// if the linear operator represents the matrix pair [A B],
// we store  pointers to A and B in this struct and then call
// sparseMat_times_mvec with A as the apply function. The apply_shifted
// function, on the other hand, will compute alpha*(A+sigma[j]*B)X+beta*Y
typedef struct TYPE(private_sparseMat_pair)
{
  TYPE(const_sparseMat_ptr) A;
  TYPE(const_sparseMat_ptr) B;
} TYPE(private_sparseMat_pair);

// same with linearOps instead of sparseMats
typedef struct TYPE(private_linearOp_pair)
{
  TYPE(const_linearOp_ptr) A;
  TYPE(const_linearOp_ptr) B;
  mutable TYPE(mvec_ptr) Xtmp;
} TYPE(private_linearOp_pair);

// same with k of linearOps
typedef struct TYPE(private_linearOp_k)
{
  TYPE(const_linearOp_ptr)* k_Ops;
  int nops;
  int * which_apply;
  const _ST_** sigma;
  int num_sigma;
  int num_shifts;
  mutable TYPE(mvec_ptr) Xtmp;
  mutable TYPE(mvec_ptr) Xtmp2;
} TYPE(private_linearOp_k);

typedef struct TYPE(private_linearOp_product)
{
  std::vector<TYPE(const_linearOp_ptr)> members_;
  mutable TYPE(mvec_ptr) Xtmp;
  mutable TYPE(mvec_ptr) Xtmp2;
} TYPE(private_linearOp_product);

//! just to have some function to point to
void SUBR(private_linearOp_destroy_nothing)(TYPE(linearOp_ptr) op, int* iflag)
{
  *iflag=0;
}

//! just to have some function to point to
void SUBR(private_linearOp_destroy_sparseMat_pair_wrapper)(TYPE(linearOp_ptr) op, int* iflag)
{
  TYPE(private_sparseMat_pair)* pair = (TYPE(private_sparseMat_pair)*)(op->A);
  delete pair; // delete the struct, not the matrices
  *iflag=0;
}

//! just to have some function to point to
void SUBR(private_linearOp_destroy_linearOp_pair_wrapper)(TYPE(linearOp_ptr) op, int* iflag)
{
  TYPE(private_linearOp_pair)* pair = (TYPE(private_linearOp_pair)*)(op->A);
  if (pair->Xtmp!=NULL) PHIST_CHK_IERR(SUBR(mvec_delete)(pair->Xtmp,iflag),*iflag);
  delete pair; // delete the struct, not the wrapped operators
  *iflag=0;
}

//! just to have some function to point to
void SUBR(private_linearOp_destroy_linearOp_k_wrapper)(TYPE(linearOp_ptr) op, int* iflag)
{
  TYPE(private_linearOp_k)* k_Ops = (TYPE(private_linearOp_k)*)(op->A);
  if (k_Ops->Xtmp!=NULL) PHIST_CHK_IERR(SUBR(mvec_delete)(k_Ops->Xtmp,iflag),*iflag);
  if (k_Ops->Xtmp2!=NULL) PHIST_CHK_IERR(SUBR(mvec_delete)(k_Ops->Xtmp2,iflag),*iflag);
  delete k_Ops; // delete the struct, not the wrapped operators
  *iflag=0;
}

// this function can be used to create an operator which encapsulates a CRS matrix.
// It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(linearOp_wrap_sparseMat)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, int* iflag)
{
  *iflag=0;
  op->A = A;
  op->aux=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&op->range_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&op->domain_map,iflag),*iflag);
  op->apply = &SUBR(sparseMat_times_mvec);
  op->applyT = &SUBR(sparseMatT_times_mvec);
  op->apply_shifted = &SUBR(sparseMat_times_mvec_vadd_mvec);
  op->fused_apply_mvTmv = &SUBR(fused_spmv_mvTmv);
  op->update=NULL;
  op->destroy=&SUBR(private_linearOp_destroy_nothing);
  return;
}

// helper function to compute Y=beta*Y + alpha*AX, ignoring B
void SUBR(private_linearOp_apply_sparseMat_pair_only_A)
(_ST_ alpha, const void* AB, TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
  // get separate pointers to A and B
  TYPE(private_sparseMat_pair) *_AB = (TYPE(private_sparseMat_pair)*)AB;
  TYPE(const_sparseMat_ptr) A = _AB->A;
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(alpha,A,X,beta,Y,iflag),*iflag);
}

// helper function to compute Y=beta*Y + alpha*(A+shifts[j]*B)X_j (apply_shifted function for matrix pair)
void SUBR(private_linearOp_apply_sparseMat_pair_shifted)
(_ST_ alpha, const void* AB, _ST_ const shifts[], TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  // get separate pointers to A and B
  TYPE(private_sparseMat_pair) *_AB = (TYPE(private_sparseMat_pair)*)AB;
  TYPE(const_sparseMat_ptr) A = _AB->A;
  TYPE(const_sparseMat_ptr) B = _AB->B;

  // start by figuring out the number of vectors (=#shifts)
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,iflag),*iflag);
  // call the fused kernel function with shift arrays shift1[], shift2[]
  _ST_ shift1[nvec], shift2[nvec];
  for (int i=0; i<nvec; i++)
  {
    shift1[i]=st::one();
    shift2[i]=shifts[i];
  }
  PHIST_CHK_IERR(SUBR(fused_spmv_pair)(alpha,shift1,A,shift2,B,X,beta,Y,iflag),*iflag);
}

// helper function to compute Y=beta*Y + alpha*A*B*X (apply function for linearOp_product)
void SUBR(private_linearOp_apply_linearOp_product)
(_ST_ alpha, const void* AB, TYPE(const_mvec_ptr) X, _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  // get separate pointers to A and B
  TYPE(private_linearOp_pair) *_AB = (TYPE(private_linearOp_pair)*)AB;
  TYPE(const_linearOp_ptr) A = _AB->A;
  TYPE(const_linearOp_ptr) B = _AB->B;

  // check if Xtmp needs reallocation
  int nvX, nvXtmp=-1;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvX,iflag),*iflag);
  if (_AB->Xtmp!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(_AB->Xtmp,&nvXtmp,iflag),*iflag);
  }

  if (nvX!=nvXtmp)
  {
    if (_AB->Xtmp!=NULL)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(_AB->Xtmp,iflag),*iflag);
      _AB->Xtmp=NULL;
    }
    PHIST_CHK_IERR(SUBR(mvec_create)(&(_AB->Xtmp),B->range_map,nvX,iflag),*iflag);
  }

  TYPE(mvec_ptr) BX = _AB->Xtmp;
  PHIST_CHK_IERR(SUBR(linearOp_apply)(st::one(),B,X,st::zero(),BX,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(linearOp_apply)(alpha,A,BX,beta,Y,iflag),*iflag);
}

//
void SUBR(linearOp_wrap_sparseMat_pair)(TYPE(linearOp_ptr) op,
        TYPE(const_sparseMat_ptr) A, TYPE(const_sparseMat_ptr) B, int *iflag)
{
  // setup maps etc.
  PHIST_CHK_IERR(SUBR(linearOp_wrap_sparseMat)(op,A,iflag),*iflag);
  TYPE(private_sparseMat_pair) *pair=new TYPE(private_sparseMat_pair);
  pair->A=A;
  pair->B=B;
  op->A=(void*)(pair);
  op->aux=NULL;
  op->apply=&SUBR(private_linearOp_apply_sparseMat_pair_only_A);
  op->apply_shifted=&SUBR(private_linearOp_apply_sparseMat_pair_shifted);
  op->update=NULL;
  op->destroy=&SUBR(private_linearOp_destroy_sparseMat_pair_wrapper);
}

void SUBR(linearOp_wrap_linearOp_product)(TYPE(linearOp_ptr) op,
TYPE(const_linearOp_ptr) A, TYPE(const_linearOp_ptr) B, int* iflag)
{
  // setup maps etc.
  TYPE(private_linearOp_pair) *pair=new TYPE(private_linearOp_pair);
  pair->A=A;
  pair->B=B;
  pair->Xtmp=NULL;
  op->A=(void*)(pair);
  op->aux=NULL;

  op->range_map=A->range_map;
  op->domain_map=B->domain_map;

  op->apply=&SUBR(private_linearOp_apply_linearOp_product);
  op->apply_shifted=NULL;
  op->update=NULL;
  op->destroy=&SUBR(private_linearOp_destroy_linearOp_pair_wrapper);
}

void SUBR(private_linearOp_product_apply)(_ST_ alpha, const void* k_op, TYPE(const_mvec_ptr) X,
 _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  TYPE(private_linearOp_product) *Op_k = (TYPE(private_linearOp_product)*)k_op;
  
  //at beginning: Xtmp=X
  if (Op_k->Xtmp!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(Op_k->Xtmp,iflag),*iflag);
    Op_k->Xtmp=NULL;
  }
  Op_k->Xtmp=(TYPE(mvec_ptr))X;

  int nvX, nvXtmp2=-1;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvX,iflag),*iflag);

  // check if Xtmp2 needs reallocation
  if (Op_k->Xtmp2!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Op_k->Xtmp2,&nvXtmp2,iflag),*iflag);
    if (nvX!=nvXtmp2)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(Op_k->Xtmp2,iflag),*iflag);
      Op_k->Xtmp2=NULL;
    }
  } 
  
  int cnt = 0;
  //iteration over operator-members
  for (auto it=(Op_k->members_).begin(); it!=(Op_k->members_).end(); it++)
  {
    TYPE(const_linearOp_ptr) op_i = *it;
    
    if(nvX!=nvXtmp2)
    {
      PHIST_CHK_IERR(SUBR(mvec_create)(&(Op_k->Xtmp2),op_i->range_map,nvX,iflag),*iflag);
    }

    PHIST_CHK_IERR(SUBR(linearOp_apply_respective)(st::one(),op_i,Op_k->Xtmp,st::zero(),Op_k->Xtmp2,iflag),*iflag);
    
    if(cnt!=0)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(Op_k->Xtmp,iflag),*iflag);
      Op_k->Xtmp=NULL;
    }
    Op_k->Xtmp = Op_k->Xtmp2;
    Op_k->Xtmp2=NULL;
    nvXtmp2=-1;
    cnt++;
  }

  // Y <- alpha*Xtmp + beta*Y
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,Op_k->Xtmp,beta,Y,iflag),*iflag);
}

void SUBR(linearOp_product_extend)(TYPE(linearOp_ptr) k_op, TYPE(const_linearOp_ptr) new_member, int* iflag)
{   
#include "phist_std_typedefs.hpp"

  TYPE(private_linearOp_product) *Op_k = (TYPE(private_linearOp_product)*)k_op->A;

  // correct range map and domain map
  if(k_op->domain_map == NULL)
  {
    if(new_member->use_transpose)
    {
      k_op->domain_map = new_member->range_map;
      k_op->range_map = new_member->domain_map;
    }
    else
    {
      k_op->domain_map = new_member->domain_map;
      k_op->range_map = new_member->range_map;  
    }
  }
  else
  {
    if(new_member->use_transpose)
    {
      k_op->range_map = new_member->domain_map;
    }
    else
    {
      k_op->range_map = new_member->range_map;
    }
  }
  
  (Op_k->members_).push_back(new_member);
}

void SUBR(private_linearOp_product_delete)(TYPE(linearOp_ptr) k_op, int* iflag)
{
  TYPE(private_linearOp_product)* Op_k = (TYPE(private_linearOp_product)*)(k_op->A);
  if (Op_k->Xtmp!=NULL) PHIST_CHK_IERR(SUBR(mvec_delete)(Op_k->Xtmp,iflag),*iflag);
  if (Op_k->Xtmp2!=NULL) PHIST_CHK_IERR(SUBR(mvec_delete)(Op_k->Xtmp2,iflag),*iflag);
  delete Op_k; // delete the struct, not the wrapped operators
  *iflag=0; 
}

void SUBR(linearOp_product_create)(TYPE(linearOp_ptr) k_op, int* iflag)
{
  TYPE(private_linearOp_product) *k_Op = new TYPE(private_linearOp_product);
  k_Op->Xtmp = NULL;
  k_Op->Xtmp2 = NULL;
  k_op->A = (void*)k_Op;
  k_op->use_transpose = NULL;
  k_op->shifts = NULL;
  k_op->apply = &SUBR(private_linearOp_product_apply);
  k_op->applyT = NULL;
  k_op->apply_shifted = NULL;
  k_op->destroy = &SUBR(private_linearOp_product_delete);
  k_op->domain_map = NULL;
  k_op->range_map = NULL;
}

//
void SUBR(private_idOp_apply)(_ST_ alpha, const void* A, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
  *iflag=0;
  PHIST_TOUCH(A)
  SUBR(mvec_add_mvec)(alpha,X,beta,Y,iflag);
}

//
void SUBR(private_idOp_fused_apply_mvTmv)(_ST_ alpha, const void* A, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, TYPE(sdMat_ptr) YtX, TYPE(sdMat_ptr) YtY, int* iflag)
{
#include "phist_std_typedefs.hpp"
  int iflag_in=*iflag;
  *iflag=0;
  PHIST_TOUCH(A)
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,X,beta,Y,iflag),*iflag);
  *iflag=iflag_in;
  if (YtX!=NULL) PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Y,X,st::zero(),YtX,iflag),*iflag);
  if (YtY!=NULL) PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Y,Y,st::zero(),YtY,iflag),*iflag);
}

//
void SUBR(private_idOp_apply_shifted)(_ST_ alpha, const void* A, _ST_ const *sigma,
        TYPE(const_mvec_ptr) X,  _ST_ beta, TYPE(mvec_ptr) Y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_TOUCH(A)
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,iflag),*iflag);
  _ST_ shifts[nvec];
  for (int i=0;i<nvec;i++) shifts[i]=(st::one()+sigma[i])*alpha;
  SUBR(mvec_vadd_mvec)(shifts,X,beta,Y,iflag);
}

// setup identity operator that returns Y=alpha*X + beta*Y
void SUBR(linearOp_identity)(TYPE(linearOp_ptr) op,
                             phist_const_map_ptr  range_map,
                             phist_const_map_ptr domain_map, int* iflag)
{
  *iflag=0;
  op->A=NULL;
  op->aux=NULL;
  op->apply = &SUBR(private_idOp_apply);
  op->applyT = &SUBR(private_idOp_apply);
  op->apply_shifted = &SUBR(private_idOp_apply_shifted);
  op->fused_apply_mvTmv = &SUBR(private_idOp_fused_apply_mvTmv);
  op->update=NULL;
  op->destroy=&SUBR(private_linearOp_destroy_nothing);
  op->range_map=range_map;
  op->domain_map=domain_map;
}

//! pointer to function for computing Y=alpha*A*X+beta*Y
void SUBR(linearOp_apply)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
  {
    PHIST_CHK_IERR(*iflag=(A_op->apply==NULL)? PHIST_BAD_CAST:0,*iflag);
    PHIST_CHK_IERR(A_op->apply(alpha,A_op->A,X,beta,Y,iflag),*iflag);
  }

//! apply transpose
 void SUBR(linearOp_applyT)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
  {
    PHIST_CHK_IERR(*iflag=(A_op->applyT==NULL)? PHIST_BAD_CAST:0,*iflag);
    PHIST_CHK_IERR(A_op->applyT(alpha,A_op->A,X,beta,Y,iflag),*iflag);
  }
 //! pointer to function for computing Y=(A-sigma[j]B)*X[j]+beta*Y[j]
 void SUBR(linearOp_apply_shifted)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op, _ST_ const * sigma,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
  {
    PHIST_CHK_IERR(*iflag=(A_op->apply_shifted==NULL)? PHIST_BAD_CAST:0,*iflag);
    PHIST_CHK_IERR(A_op->apply_shifted(alpha,A_op->A,sigma,X,beta,Y,iflag),*iflag);
  }
 //! pointer to function for computing Y=alpha*A*X+beta*Y
 void SUBR(linearOp_apply_respective)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op,
        TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* iflag)
  {
    if(A_op->shifts != NULL)
    {  
      if(A_op->use_transpose)
      {
        PHIST_SOUT(PHIST_ERROR, "not clear whether applyT or apply_shifted should be used \n",
                              "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
        *iflag=-1;
        return;
      }
      else
      {
        PHIST_SOUT(PHIST_VERBOSE, "linearOp_apply_respective uses linearOp_apply_shifted \n",
                              "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
        PHIST_CHK_IERR(SUBR(linearOp_apply_shifted)(alpha,A_op,A_op->shifts,X,beta,Y,iflag),*iflag);
      }
    }  
    
    else if(A_op->use_transpose)
    {
        PHIST_SOUT(PHIST_VERBOSE, "linearOp_apply_respective uses linearOp_applyT \n",
                              "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
        PHIST_CHK_IERR(SUBR(linearOp_applyT)(alpha,A_op,X,beta,Y,iflag),*iflag);
    }
    
    else
    {
        PHIST_SOUT(PHIST_VERBOSE, "linearOp_apply_respective uses linearOp_apply \n",
                              "(file %s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
        PHIST_CHK_IERR(SUBR(linearOp_apply)(alpha,A_op,X,beta,Y,iflag),*iflag);
    }
  }
//! apply operator and compute inner products with in- and output vector
  void SUBR(linearOp_fused_apply_mvTmv)(_ST_ alpha, TYPE(const_linearOp_ptr) A_op, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                 TYPE(mvec_ptr)        W,
                            TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                            int* iflag)
  {
    PHIST_CHK_IERR(*iflag=(A_op->fused_apply_mvTmv==NULL)? PHIST_BAD_CAST:0,*iflag);
    PHIST_CHK_IERR(A_op->fused_apply_mvTmv(alpha,A_op->A,V,beta,W,WtW,VtW,iflag),*iflag);
  }


  void SUBR(linearOp_update)(TYPE(linearOp_ptr) A_op, _ST_ sigma,
                        TYPE(const_mvec_ptr) Vkern,
                        TYPE(const_mvec_ptr) BVkern,
                        int *iflag)
  {
    PHIST_CHK_IERR(A_op->update(A_op->A,A_op->aux,sigma,Vkern,BVkern,iflag),*iflag);
  }

  //! this function can be used to clean up any data the operator may *own*,
  //! if the operator is just a wrapper for some other object that is created
  //! and deleted separately, this function should not do anything.
  //! The me object itself should *not* be free'd.
  void SUBR(linearOp_destroy)(TYPE(linearOp_ptr) A_op, int* iflag)
  {
    if (A_op->destroy!=NULL)
    {
      PHIST_CHK_IERR(A_op->destroy(A_op,iflag),*iflag);
    }
  }

} // extern "C"

