extern "C" {

// this function can be used to create an operator which encapsulates a CRS matrix.
// It does not allocate memory for the op struct, the caller has to do that beforehand.
void SUBR(op_wrap_crsMat)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A, int* ierr)
  {
  *ierr=0;
  op->A = A;
  PHIST_CHK_IERR(SUBR(crsMat_get_range_map)(A,&op->range_map,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(crsMat_get_domain_map)(A,&op->domain_map,ierr),*ierr);
  op->apply = &SUBR(crsMat_times_mvec);
  op->applyT = &SUBR(crsMatT_times_mvec);
  op->apply_shifted = &SUBR(crsMat_times_mvec_vadd_mvec);
  return;
  }


//
void SUBR(private_idOp_apply)(_ST_ alpha, const void* A, TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
{
  *ierr=0;
  TOUCH(A)
  SUBR(mvec_add_mvec)(alpha,X,beta,Y,ierr);
}

//
void SUBR(private_idOp_apply_shifted)(_ST_ alpha, const void* A, _ST_ const *sigma, 
        TYPE(const_mvec_ptr) X,  _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
{
  *ierr=0;
  TOUCH(A)
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,ierr),*ierr);
  _ST_ shifts[nvec];
  for (int i=0;i<nvec;i++) shifts[i]=sigma[i]*alpha;
  SUBR(mvec_vadd_mvec)(shifts,X,beta,Y,ierr);
}

// setup identity operator that returns Y=alpha*X + beta*Y
void SUBR(op_identity)(TYPE(op_ptr) op, int* ierr)
{
  *ierr=0;
  op->A=NULL;
  op->apply = &SUBR(private_idOp_apply);
  op->applyT = &SUBR(private_idOp_apply);
  op->apply_shifted = &SUBR(private_idOp_apply_shifted);
}

typedef struct TYPE(private_carpOpData)
{
  TYPE(carpData)* carp_;
  TYPE(const_crsMat_ptr) A_;
  TYPE(const_mvec_ptr) B_;
} TYPE(private_carpOpData);


//
void SUBR(private_I_minus_dkswp_shifted)(_ST_ alpha, const void* vdat, _ST_ const *shifts,
        TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
{
#include "phist_std_typedefs.hpp"
ENTER_FCN(__FUNCTION__);
CAST_PTR_FROM_VOID(TYPE(private_carpOpData) const, dat, vdat,*ierr);

  if (alpha!=st::one() || beta!=st::zero())
  {
    PHIST_SOUT(PHIST_ERROR,"only case alpha=1, beta=0 implemented in %s\n"
                           "(%s, line %d)\n",__FUNCTION__,__FILE__,__LINE__);
    *ierr=-99;
    return;
  }
  // TODO - optimize this implementation of I-DKSWP(A)
//  TYPE(mvec_ptr) X_tmp=NULL;
//  PHIST_CHK_IERR(SUBR(mvec_clone)(&X_tmp,X,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(carp_fb)(dat->carp_, dat->A_, dat->B_, shifts,
                X, Y, ierr), *ierr);

//  SUBR(mvec_add_mvec)(alpha,X,beta,Y,ierr);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),X,-st::one(),Y,ierr),*ierr);
}

void SUBR(private_I_minus_dkswp)(_ST_ alpha, const void* vdat,
        TYPE(const_mvec_ptr) X,
        _ST_ beta, TYPE(mvec_ptr) Y, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nvec,ierr),*ierr);
  _ST_ shifts[nvec];
  for (int i=0;i<nvec;i++) shifts[i]=st::one();
  PHIST_CHK_IERR(SUBR(private_I_minus_dkswp_shifted)(alpha,vdat,shifts,X,beta,Y,ierr),*ierr);
}

//! given a crsMat, create the operator I-DCSWP(A,x), double CARP sweep.
//! This operator can be passed to CG for solving Ax=b via the normal
//! equations AA'y=b, x=A'y. We allow for an array of shifts here, so
//! that the operator acting on each column i is in fact
//! I-DCSWP(A-s[i]I,x[i]). If shift is NULL, s=0 is assumed.
void SUBR(op_carp)(TYPE(op_ptr) op, TYPE(const_crsMat_ptr) A,
        _MT_ omega, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  // note: this operator acts as A*A^T
  PHIST_CHK_IERR(SUBR(crsMat_get_range_map)(A,&op->range_map,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(crsMat_get_range_map)(A,&op->domain_map,ierr),*ierr);
  
  TYPE(private_carpOpData)* dat =
        (TYPE(private_carpOpData)*)malloc(sizeof(TYPE(private_carpOpData)));
  dat->A_=A;
  dat->B_=NULL;
  PHIST_CHK_IERR(SUBR(carp_create)(A,&(dat->carp_),ierr),*ierr);
  dat->carp_->omega_=omega;
  op->A=dat;
  op->apply=&SUBR(private_I_minus_dkswp);
  op->applyT=NULL; // not supported
  op->apply_shifted = &SUBR(private_I_minus_dkswp_shifted);

  return;
}

} // extern "C"
