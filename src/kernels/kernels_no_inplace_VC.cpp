//! if a kernel lib does not provide a specialized implementation
//! for the operation V <- V*C it can use this default one by including
//! this source file
extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, 
TYPE(const_sdMat_ptr) M_, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;

  lidx_t chunkSize = 64;

  // get dimensions
  lidx_t nV;
  int nvec, nM, mM;
  PHIST_CHK_IERR( SUBR( mvec_my_length     ) (V_, &nV,        iflag), *iflag);
  PHIST_CHK_IERR( SUBR( mvec_num_vectors   ) (V_, &nvec,      iflag), *iflag);
  PHIST_CHK_IERR( SUBR( sdMat_get_nrows    ) (M_, &nM,        iflag), *iflag);
  PHIST_CHK_IERR( SUBR( sdMat_get_ncols    ) (M_, &mM,        iflag), *iflag);

  // check dimensions
  PHIST_CHK_IERR(*iflag = (nvec == nM ? 0 : -1), *iflag);
  PHIST_CHK_IERR(*iflag = (nvec >= mM ? 0 : -1), *iflag);

  // get map
  const_map_ptr_t map = NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V_, &map, iflag), *iflag);

  // create temporary mvec
  TYPE(mvec_ptr) Vtmp = NULL;
  PHIST_CHK_IERR(SUBR(mvec_create)(&Vtmp, map, mM, iflag), *iflag);

  // calculate and set block
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(), V_, M_, st::zero(), Vtmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_set_block)(V_, Vtmp, 0, mM-1, iflag), *iflag);

  // delete temp data
  PHIST_CHK_IERR(SUBR(mvec_delete)(Vtmp, iflag), *iflag);
}

