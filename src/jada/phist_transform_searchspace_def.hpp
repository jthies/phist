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

