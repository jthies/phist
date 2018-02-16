/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
extern "C" void SUBR(fused_spmv_pair)(_ST_ alpha, 
                                       _ST_ const shift1[], TYPE(const_sparseMat_ptr) A1, 
                                       _ST_ const shift2[], TYPE(const_sparseMat_ptr) A2,
                                       TYPE(const_mvec_ptr) X,
                            _ST_ beta, TYPE(mvec_ptr)       Y,
                            int* iflag)
                                                                                    
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_WARN_MISSING_KERNEL(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  // pass *iflag to both kernels.
  int iflag0=0, iflag1=0, iflag2=0, iflag3=0, iflag4=0,iflag5=0, iflag6=0, iflag7=0;
  // we need a temporary vector because there is no function for shift[j]*A*X_j
  int nvec=-1;
  SUBR(mvec_num_vectors)(X,&nvec,&iflag0);
  phist_const_map_ptr map=NULL;
  SUBR(mvec_get_map)(X,&map,&iflag1);
  TYPE(mvec_ptr) tmp=NULL;
  SUBR(mvec_create)(&tmp,map,nvec,&iflag2);
  iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A1,X,st::zero(),tmp,&iflag3);
  SUBR(mvec_vadd_mvec)(shift1,tmp,beta,Y,&iflag4);

  iflag1 = *iflag;
  SUBR(sparseMat_times_mvec)(alpha,A2,X,st::zero(),tmp,&iflag5);
  SUBR(mvec_vadd_mvec)(shift2,tmp,st::one(),Y,&iflag6);
  

  SUBR(mvec_delete)(tmp,&iflag7);

  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  PHIST_CHK_IERR(*iflag = iflag3,*iflag);
  PHIST_CHK_IERR(*iflag = iflag4,*iflag);
  PHIST_CHK_IERR(*iflag = iflag5,*iflag);
  PHIST_CHK_IERR(*iflag = iflag6,*iflag);
  PHIST_CHK_IERR(*iflag = iflag7,*iflag);
}


