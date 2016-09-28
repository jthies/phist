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

// this function can be used to create an operator which encapsulates a CRS matrix.
// It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(linearOp_wrap_sparseMat)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, int* iflag)
{
  *iflag=0;
  op->A = A;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&op->range_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&op->domain_map,iflag),*iflag);
  op->apply = &SUBR(sparseMat_times_mvec);
  op->applyT = &SUBR(sparseMatT_times_mvec);
  op->apply_shifted = &SUBR(sparseMat_times_mvec_vadd_mvec);
  op->fused_apply_mvTmv = &SUBR(fused_spmv_mvTmv);
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
  op->apply=&SUBR(private_linearOp_apply_sparseMat_pair_only_A);
  op->apply_shifted=&SUBR(private_linearOp_apply_sparseMat_pair_shifted);
  op->destroy=&SUBR(private_linearOp_destroy_sparseMat_pair_wrapper);
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
  *iflag=0;
  PHIST_TOUCH(A)
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,iflag),*iflag);
  _ST_ shifts[nvec];
  for (int i=0;i<nvec;i++) shifts[i]=sigma[i]*alpha;
  SUBR(mvec_vadd_mvec)(shifts,X,beta,Y,iflag);
}

// setup identity operator that returns Y=alpha*X + beta*Y
void SUBR(linearOp_identity)(TYPE(linearOp_ptr) op, int* iflag)
{
  *iflag=0;
  op->A=NULL;
  op->apply = &SUBR(private_idOp_apply);
  op->applyT = &SUBR(private_idOp_apply);
  op->apply_shifted = &SUBR(private_idOp_apply_shifted);
  op->fused_apply_mvTmv = &SUBR(private_idOp_fused_apply_mvTmv);
  op->destroy=&SUBR(private_linearOp_destroy_nothing);
}

} // extern "C"
