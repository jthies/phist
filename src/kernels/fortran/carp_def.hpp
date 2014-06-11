#include "phist_carp_decl.h"

extern "C" {

void SUBR(carp_setup)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ **nrms_ai2i, void** work, int* ierr)
{
  *ierr=-99;
  return;
}


void SUBR(dkswp)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs, 
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* ierr)
{
  *ierr=-99;
  return;
}

void SUBR(carp_destroy)(TYPE(const_crsMat_ptr) A, int numShifts,
_MT_* nrms_ai2i, void* work, int *ierr)
{
  *ierr=-99;
  return;
}






} // extern "C"
