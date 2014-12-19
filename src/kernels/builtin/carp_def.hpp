#include "phist_carp_decl.h"

extern "C" {

void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ **nrms_ai2i, void** work, int* iflag)
{
  void SUBR(carp_setup_f)(TYPE(const_sparseMat_ptr) A, int numShifts,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ *nrms_ai2i, void** work, int* iflag);

  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  const_map_ptr_t map;
  int nlocal;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nlocal,iflag),*iflag);
  *nrms_ai2i=new _MT_[nlocal*numShifts]; // column-major array with
                                         // nrms_ai2i(i,j) the inverse
                                         // 2-norm of row i for shift j
  PHIST_CHK_IERR(SUBR(carp_setup_f)
        (A,numShifts,sigma_r,sigma_i,*nrms_ai2i,work,iflag),*iflag);
  return;
}


void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* iflag)
{
  void SUBR(carp_sweep_f)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* iflag);
  void SUBR(carp_sweep_real_f)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* iflag);
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  bool is_complex=!(X_i==NULL||sigma_i==NULL);
  if (is_complex)
  {
    is_complex=false;
    for (int i=0;i<numShifts;i++)
    {
      if (sigma_i[i]!=mt::zero()) 
      {
        is_complex=true;
      }
    }
  }
  if (is_complex)
  {
    PHIST_CHK_IERR(SUBR(carp_sweep_f)(A, numShifts,
        sigma_r, sigma_i, Rhs, X_r, X_i, nrm_ai2i, work,
        omega, iflag),*iflag);
  }
  else
  {
    PHIST_CHK_IERR(SUBR(carp_sweep_real_f)(A, numShifts,
        sigma_r, Rhs, X_r, nrm_ai2i, work,
        omega, iflag),*iflag);
  }
  return;
}

void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A, int numShifts,
_MT_* nrms_ai2i, void* work, int *iflag)
{
  void SUBR(carp_destroy_f)(TYPE(const_sparseMat_ptr) A, int numShifts,
        void* work, int *iflag);

  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  delete [] nrms_ai2i;
  PHIST_CHK_IERR(SUBR(carp_destroy_f)(A,numShifts,work,iflag),*iflag);
  return;
}


} // extern "C"
