void SUBR(fused_spmv_mvTmv)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                             _ST_ beta,                               TYPE(mvec_ptr)        W,
                             TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW, TYPE(sdMat_ptr) W0tW,
                             int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call separate kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag,iflag2=0,iflag3=0,iflag4=0,iflag5=0;
  TYPE(mvec_ptr) W0=NULL;
  if (W0tW!=NULL)
  {
    phist_const_map_ptr map=NULL;
      int nvW;
    SUBR(mvec_get_map)(W,&map,&iflag4);
    if (!iflag4) SUBR(mvec_num_vectors)(W,&nvW,&iflag4);
    if (!iflag4) SUBR(mvec_create)(&W0,map,nvW,&iflag4);
    if (!iflag4) SUBR(mvec_add_mvec)(st::one(),W,st::zero(),W0,&iflag4);
  }
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  if (iflag4==0 && W0tW!=NULL)
  {
    iflag4 = *iflag;
    SUBR(mvecT_times_mvec)(st::one(),W0,W,st::zero(),W0tW,&iflag4);
    SUBR(mvec_delete)(W0,&iflag5);
  }
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
  PHIST_CHK_IERR(*iflag = iflag4,*iflag);
  PHIST_CHK_IERR(*iflag = iflag5,*iflag);
}
