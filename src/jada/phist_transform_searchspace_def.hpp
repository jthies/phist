//! apply a transformation matrix M to a given search space V and AV, BV and the projection H=V'AV
void SUBR(transform_searchSpace)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) M, bool generalizedEigenproblem, int *ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // get dimensions to create temporary storage
  int nV, minBase;
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(SUBR( sdMat_get_nrows ) (M, &nV, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_get_ncols ) (M, &minBase, ierr), *ierr);

  // update V, AV and BV
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (V,  M, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (AV, M, ierr), *ierr);
  if( generalizedEigenproblem )
  {
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (BV, M, ierr), *ierr);
  }

  // we need some communicator for ghost...
  PHIST_CHK_IERR(SUBR( mvec_get_comm ) (V, &comm, ierr), *ierr);
  TYPE(sdMat_ptr) H_ = NULL;
  TYPE(sdMat_ptr) Htmp = NULL;
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp, nV, minBase, comm, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H,    &H_,    0, minBase-1,    0, minBase-1,     ierr), *ierr);

  // update H <- M' * H * M
  PHIST_CHK_IERR(SUBR( sdMat_times_sdMat )(st::one(), H, M,    st::zero(), Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat)(st::one(), M, Htmp, st::zero(), H_,   ierr), *ierr);

/*
  // set M <- I
  _ST_ *Mraw = NULL;
  lidx_t ldaM;
  PHIST_CHK_IERR(SUBR( sdMat_extract_view )(M, &Mraw, &ldaM, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_put_value )(M, st::zero(), ierr), *ierr);
  for(int i = 0; i < minBase; i++)
    Mraw[ldaM*i+i] = st::one();
*/
  // delete temp. storage
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,   ierr), *ierr);
}


#ifndef PHIST_KERNEL_LIB_FORTRAN
//! in order to shrink a subspace we need a fast in-place operation for mvec_times_sdMat
void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, TYPE(sdMat_ptr) M_, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
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
#endif
