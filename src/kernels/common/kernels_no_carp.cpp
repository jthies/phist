/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

extern "C" void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _ST_ const sigma[],
        void** work, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}


extern "C" void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A, 
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

extern "C" void SUBR(carp_sweep_aug)(TYPE(const_sparseMat_ptr) A,
        _ST_ const sigma_r[],
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

#ifndef IS_COMPLEX

extern "C" void SUBR(carp_setup_rc)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}


extern "C" void SUBR(carp_sweep_rc)(TYPE(const_sparseMat_ptr) A, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

extern "C" void SUBR(carp_sweep_aug_rc)(TYPE(const_sparseMat_ptr) A,
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

#endif


extern "C" void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A,
        void* work, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}



