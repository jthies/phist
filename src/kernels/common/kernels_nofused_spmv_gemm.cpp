void SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_self)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                            _ST_ beta,                         TYPE(mvec_ptr)        W,
                                                                                               TYPE(sdMat_ptr)       D,
                                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}


void SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_other)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                            _ST_ beta,                         TYPE(mvec_ptr)        W,
                                                                                               TYPE(sdMat_ptr)       D,
                                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,V,st::zero(),D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}


void SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_both)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                            _ST_ beta,                         TYPE(mvec_ptr)        W,
                                                                        TYPE(sdMat_ptr)     C, TYPE(sdMat_ptr)       D,
                                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,V,st::zero(),C,&iflag2);
  int iflag3 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),D,&iflag3);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  PHIST_CHK_IERR(*iflag = iflag3,*iflag);
}

