void SUBR(mvec_times_sdMat_augmented)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                  TYPE(const_sdMat_ptr) C,
                                      _ST_ beta,  TYPE(mvec_ptr)        W,
                                                  TYPE(sdMat_ptr)       D,
                                                  int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // remember flags
  int flags = *iflag;
  // call two kernels
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(alpha,V,C,beta,W,iflag),*iflag);
  *iflag = flags;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),D,iflag),*iflag);
}

void SUBR(mvecT_times_mvec_times_sdMat_inplace)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                            TYPE(mvec_ptr)        W,
                                                            TYPE(const_sdMat_ptr) C,
                                                _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // remember flags
  int flags = *iflag;
  // call two kernels
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(W,C,iflag),*iflag);
  *iflag = flags;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(alpha,V,W,beta,D,iflag),*iflag);
}

void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V, 
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


