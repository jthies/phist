//! apply a transformation matrix M to a given search space V and AV, BV and the projection H=V'AV
void SUBR(transform_searchSpace)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) M, bool transformBV, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;

  // get dimensions to create temporary storage
  int nV, minBase;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(SUBR( sdMat_get_nrows ) (M, &nV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_get_ncols ) (M, &minBase, iflag), *iflag);

  // update V, AV and BV
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (V,  M, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (AV, M, iflag), *iflag);
  if( transformBV )
  {
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (BV, M, iflag), *iflag);
  }

  // we need some communicator for ghost...
  PHIST_CHK_IERR(SUBR( mvec_get_comm ) (V, &comm, iflag), *iflag);
  TYPE(sdMat_ptr) H_ = NULL;
  TYPE(sdMat_ptr) Htmp = NULL;
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp, nV, minBase, comm, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H,    &H_,    0, minBase-1,    0, minBase-1,     iflag), *iflag);

  // update H <- M' * H * M
  PHIST_CHK_IERR(SUBR( sdMat_times_sdMat )(st::one(), H, M,    st::zero(), Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat)(st::one(), M, Htmp, st::zero(), H_,   iflag), *iflag);

/*
  // set M <- I
  _ST_ *Mraw = NULL;
  phist_lidx ldaM;
  PHIST_CHK_IERR(SUBR( sdMat_extract_view )(M, &Mraw, &ldaM, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_put_value )(M, st::zero(), iflag), *iflag);
  for(int i = 0; i < minBase; i++)
    Mraw[ldaM*i+i] = st::one();
*/
  // delete temp. storage
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,   iflag), *iflag);
}

//! apply transformation matrices MV and MW to given search spaces V and W, resp.
//! Update BV and the projections H=W'V, H_A=W'AV
void SUBR(transform_searchSpaceHarmonic)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) W, TYPE(mvec_ptr) BV,
        TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) H_A, TYPE(sdMat_ptr) MV, TYPE(sdMat_ptr) MW, bool transformBV, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,MV,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(W,MW,iflag),*iflag);
  if (transformBV)
  {
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(BV,MV,iflag),*iflag);
  }

  // we need some communicator for ghost...
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(SUBR( mvec_get_comm ) (V, &comm, iflag), *iflag);

  int nV, minBase;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V, &nV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(MV, &minBase, iflag), *iflag);
  
  TYPE(sdMat_ptr) H_ = NULL;
  TYPE(sdMat_ptr) Htmp = NULL;
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp, nV, minBase, comm, iflag), *iflag);

  // update H <- MV' * H * MW
  PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H,    &H_,    0, minBase-1,    0, minBase-1,     iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_times_sdMat )(st::one(), H, MW,    st::zero(), Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat)(st::one(), MV, Htmp, st::zero(), H_,   iflag), *iflag);

  // update H_A <- MV' * H_A * MW
  PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H_A,    &H_,    0, minBase-1,    0, minBase-1,     iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_times_sdMat )(st::one(), H_A, MW,    st::zero(), Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat)(st::one(), MV, Htmp, st::zero(), H_,   iflag), *iflag);

  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp, iflag), *iflag);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,   iflag), *iflag);

  return;
}
