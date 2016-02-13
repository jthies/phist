extern "C" void SUBR(sparseMat_times_mvec_fused_norm2)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                            _ST_ beta,                         TYPE(mvec_ptr)        W,
                                                                               _MT_*                 Wnrm,
                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvec_norm2)(W,Wnrm,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}


extern "C" void SUBR(sparseMat_times_mvec_fused_dot)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                          _ST_ beta,                         TYPE(mvec_ptr)        W,
                                                                             _ST_*                 WdotV,
                                          int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvec_dot_mvec)(W,V,WdotV,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}


extern "C" void SUBR(sparseMat_times_mvec_fused_dot_norm2)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                _ST_ beta,                         TYPE(mvec_ptr)        W,
                                                            _ST_*           WdotV, _MT_*                 Wnrm,
                                                int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvec_dot_mvec)(W,V,WdotV,&iflag2);
  int iflag3 = *iflag;
  SUBR(mvec_norm2)(W,Wnrm,&iflag3);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  PHIST_CHK_IERR(*iflag = iflag3,*iflag);
}

