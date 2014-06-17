#include "phist_carp_decl.h"

extern "C" {

void SUBR(carp_setup)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ **nrms_ai2i, void** work, int* ierr)
{
  void SUBR(carp_setup_f)(TYPE(const_crsMat_ptr) A, int numShifts,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ **nrms_ai2i, void** work, int* ierr);

  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  PHIST_CHK_IERR(SUBR(carp_setup_f)
        (A,numShifts,sigma_r,sigma_i,nrms_ai2i,work,ierr),*ierr);
  return;
}


void SUBR(carp_sweep)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* ierr)
{
  void SUBR(carp_sweep_f)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* ierr);
  
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  PHIST_CHK_IERR(SUBR(carp_sweep_f)(A, numShifts,
        sigma_r, sigma_i, Rhs, X_r, X_i, nrm_ai2i, work,
        omega, ierr),*ierr);
  return;
}

void SUBR(carp_destroy)(TYPE(const_crsMat_ptr) A, int numShifts,
_MT_* nrms_ai2i, void* work, int *ierr)
{
  void SUBR(carp_destroy_f)(TYPE(const_crsMat_ptr) A, int numShifts,
        _MT_* nrms_ai2i, void* work, int *ierr);

  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  PHIST_CHK_IERR(SUBR(carp_destroy_f)(A,numShifts,nrms_ai2i,work,ierr),*ierr);
  return;
}


} // extern "C"
