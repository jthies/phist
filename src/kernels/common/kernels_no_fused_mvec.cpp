/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
extern "C" void SUBR(fused_mvsd_mvTmv)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                  TYPE(const_sdMat_ptr) C,
                                      _ST_ beta,  TYPE(mvec_ptr)        W,
                                                  TYPE(sdMat_ptr)       WtW,
                                                  int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  // pass *iflag to both kernels, it may e.g. contain PHIST_NO_GLOBAL_REDUCTION
  int iflag1 = *iflag;
  SUBR(mvec_times_sdMat)(alpha,V,C,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),WtW,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

extern "C" void SUBR(fused_mvsdi_mvTmv)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                            TYPE(mvec_ptr)        W,
                                                            TYPE(const_sdMat_ptr) C,
                                                _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  // pass *iflag to both kernels, it may e.g. contain PHIST_NO_GLOBAL_REDUCTION
  int iflag1 = *iflag;
  SUBR(mvec_times_sdMat_inplace)(W,C,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(alpha,V,W,beta,D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

