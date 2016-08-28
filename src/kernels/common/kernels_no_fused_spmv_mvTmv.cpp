void SUBR(fused_spmv_mvTmv)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                             _ST_ beta,                               TYPE(mvec_ptr)        W,
                             TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                             int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  TYPE(mvec_ptr) _W=W;
  _ST_ _alpha=alpha, _beta=beta;
  int iflag_in=*iflag;
  *iflag=0;
  
  // pathological case: no output args at all
  if (WtW==NULL && VtW==NULL && W==NULL) return;
  // use temporary vector if necessary (this can of course be avoided in an optimized implementation)
  if (W==NULL && (WtW!=NULL || VtW!=NULL))
  {
    phist_const_map_ptr map=NULL;
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(W,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nvec,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&_W,map,nvec,iflag),*iflag);
    _alpha=st::one();
    _beta=st::zero();
  }
  
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = iflag_in,iflag2=0,iflag3=0;
  SUBR(sparseMat_times_mvec)(_alpha,A,V,_beta,_W,&iflag1);
  if (WtW!=NULL)
  {
    iflag2 = iflag_in;
    _ST_ scal =st::one()/(_alpha*_alpha);
    SUBR(mvecT_times_mvec)(st::one(),_W,_W,st::zero(),WtW,&iflag2);
  }
  if (VtW!=NULL)
  {
    iflag3 = iflag_in;
    _ST_ scal =st::one()/(_alpha);
    SUBR(mvecT_times_mvec)(st::one(),V,_W,st::zero(),VtW,&iflag2);
  }
  if (W==NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(_W,iflag),*iflag);
  }
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  PHIST_CHK_IERR(*iflag = iflag3,*iflag);
}
