//! if a kernel lib does not provide a specialized implementation
//! for the operation V <- V*C it can use this default one by including
//! this source file
extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, 
TYPE(const_sdMat_ptr) M_, int* ierr)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  lidx_t chunkSize = 64;

  // get dimensions
  lidx_t nV;
  int nvec, nM, mM;
  PHIST_CHK_IERR( SUBR( mvec_my_length     ) (V_, &nV,        ierr), *ierr);
  PHIST_CHK_IERR( SUBR( mvec_num_vectors   ) (V_, &nvec,      ierr), *ierr);
  PHIST_CHK_IERR( SUBR( sdMat_get_nrows    ) (M_, &nM,        ierr), *ierr);
  PHIST_CHK_IERR( SUBR( sdMat_get_ncols    ) (M_, &mM,        ierr), *ierr);

  // check dimensions
  PHIST_CHK_IERR(*ierr = (nvec == nM ? 0 : -1), *ierr);
  PHIST_CHK_IERR(*ierr = (nvec >= mM ? 0 : -1), *ierr);

  // get map
  const_map_ptr_t map = NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V_, &map, ierr), *ierr);

  // create temporary mvec
  TYPE(mvec_ptr) Vtmp = NULL;
  PHIST_CHK_IERR(SUBR(mvec_create)(&Vtmp, map, mM, ierr), *ierr);

  // calculate and set block
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(), V_, M_, st::zero(), Vtmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V_, Vtmp, 0, mM-1, ierr), *ierr);

  // delete temp data
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vtmp, ierr), *ierr);
}

