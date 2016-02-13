extern "C" void SUBR(mvec_times_sdMat_aug)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                  TYPE(const_sdMat_ptr) C,
                                      _ST_ beta,  TYPE(mvec_ptr)        W,
                                                  TYPE(sdMat_ptr)       D,
                                                  int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(mvec_times_sdMat)(alpha,V,C,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

extern "C" void SUBR(mvecT_times_mvec_times_sdMat_inplace)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                            TYPE(mvec_ptr)        W,
                                                            TYPE(const_sdMat_ptr) C,
                                                _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(mvec_times_sdMat_inplace)(W,C,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(alpha,V,W,beta,D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

extern "C" void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V, 
                                                 TYPE(const_sdMat_ptr) C,
                                                 TYPE(mvec_ptr) W, 
                                                 TYPE(const_sdMat_ptr) D,
                                                 int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // remember flags
  int flags = *iflag;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(W,D,iflag),*iflag);
  *iflag = flags;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),V,C,st::one(),W,iflag),*iflag);
}


