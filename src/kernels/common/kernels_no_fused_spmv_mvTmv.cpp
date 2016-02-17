void SUBR(fused_spmv_mvTmv)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                             _ST_ beta,                               TYPE(mvec_ptr)        W,
                             TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                             int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag,iflag2=0,iflag3=0;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  if (WtW!=NULL)
  {
    iflag2 = *iflag;
    SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),WtW,&iflag2);
  }
  if (VtW!=NULL)
  {
    iflag3 = *iflag;
    SUBR(mvecT_times_mvec)(st::one(),V,W,st::zero(),VtW,&iflag2);
  }
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  PHIST_CHK_IERR(*iflag = iflag3,*iflag);
}
