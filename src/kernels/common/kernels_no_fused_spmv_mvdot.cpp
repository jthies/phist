/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
extern "C" void SUBR(fused_spmv_mvdot)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                                     TYPE(mvec_ptr)  W,
                            _ST_* WdotW, _ST_* VdotW,
                            int* iflag)
                                                                                    
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_WARN_MISSING_KERNEL(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  // pass *iflag to both kernels.
  int iflag1=0, iflag2=0, iflag3=0;
  iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A,V,beta,W,&iflag1);
  if (VdotW!=NULL)
  {
    iflag2 = *iflag;
    SUBR(mvec_dot_mvec)(V,W,VdotW,&iflag2);
  }
  if (WdotW!=NULL)
  {
    iflag3 = *iflag;
    SUBR(mvec_dot_mvec)(W,W,WdotW,&iflag3);
  }
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  PHIST_CHK_IERR(*iflag = iflag3,*iflag);
}


extern "C" void SUBR(fused_spmv_mvdot_mvadd)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                                     TYPE(mvec_ptr)  W,
                            _ST_ gamma, _ST_ delta,                        TYPE(mvec_ptr)  U,
                            _ST_* WdotW, _ST_* VdotW,
                            int* iflag)
                                                                                    
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_WARN_MISSING_KERNEL(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  // pass *iflag to both kernels.
  int iflag1=0, iflag2=0;
  
  iflag1=*iflag;
  SUBR(fused_spmv_mvdot)(alpha,A,V,beta,W,WdotW,VdotW,&iflag1);
  if (U!=NULL)
  {
    iflag2=*iflag;
    SUBR(mvec_add_mvec)(gamma,W,delta,U,&iflag2);
  }

  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

