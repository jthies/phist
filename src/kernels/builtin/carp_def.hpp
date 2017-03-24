/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

extern "C" {

void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _ST_ const sigma[],
        void** work, int* iflag)
{
  void SUBR(carp_setup_f)(TYPE(const_sparseMat_ptr) A, int numShifts,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag);

  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  phist_const_map_ptr map;
  int nlocal;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(carp_setup_f)
        (A,numShifts,sigma,NULL,work,iflag),*iflag);
  return;
}


void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A, _ST_ const sigma[],
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X,
        void* const work, _MT_ const * omega, int* iflag)
{
  void SUBR(carp_sweep_real_f)(TYPE(const_sparseMat_ptr) A, 
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r,
        void* const work,
        _MT_ const * omega, int* iflag);
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CHK_IERR(SUBR(carp_sweep_real_f)(A,
      sigma, Rhs, X, work,
      omega, iflag),*iflag);
  return;
}

void SUBR(carp_setup_rc)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
  void SUBR(carp_setup_f)(TYPE(const_sparseMat_ptr) A, int numShifts,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag);

  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  phist_const_map_ptr map;
  int nlocal;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(carp_setup_f)
        (A,numShifts,sigma_r,sigma_i,work,iflag),*iflag);
  return;
}


void SUBR(carp_sweep_rc)(TYPE(const_sparseMat_ptr) A, _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work, _MT_ const * omega, int* iflag)
{
  void SUBR(carp_sweep_f)(TYPE(const_sparseMat_ptr) A, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag);

#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CHK_IERR(SUBR(carp_sweep_f)(A,
      sigma_r, sigma_i, Rhs, X_r, X_i, work,
      omega, iflag),*iflag);
  return;
}

void SUBR(carp_sweep_aug)(TYPE(const_sparseMat_ptr) A,
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) Q,
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X,
        TYPE(sdMat_ptr) q,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

void SUBR(carp_sweep_aug_rc)(TYPE(const_sparseMat_ptr) A,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Q,
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        TYPE(sdMat_ptr) q_r, TYPE(sdMat_ptr) q_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A, void* work, int *iflag)
{
  void SUBR(carp_destroy_f)(TYPE(const_sparseMat_ptr) A, void* work, int *iflag);

  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CHK_IERR(SUBR(carp_destroy_f)(A,work,iflag),*iflag);
  return;
}


} // extern "C"
