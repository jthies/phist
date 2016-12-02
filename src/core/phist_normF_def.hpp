void SUBR(sdMat_normF)(TYPE(const_sdMat_ptr) M, _MT_ *f, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  int n = 0;
  TYPE(sdMat_ptr) tmp = NULL;
  _ST_ *raw_tmp = NULL;
  phist_lidx lda_tmp = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp,n,n,NULL,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(),M,M,st::zero(),tmp,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(tmp,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp,&raw_tmp,&lda_tmp,iflag),*iflag);
  *f = mt::zero();
  for(int i = 0; i < n; i++)
    *f += st::real(raw_tmp[i*lda_tmp+i]);
  *f = mt::sqrt(*f);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp,iflag),*iflag);
}

void SUBR(mvec_normF)(TYPE(const_mvec_ptr) V, _MT_ *f, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  int n = 0;
  int flags = *iflag;
  TYPE(sdMat_ptr) tmp = NULL;
  _ST_ *raw_tmp = NULL;
  phist_lidx lda_tmp = 0;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&tmp,n,n,NULL,iflag),*iflag);
  *iflag = flags;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),tmp,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(tmp,&raw_tmp,&lda_tmp,iflag),*iflag);
  *f = mt::zero();
  for(int i = 0; i < n; i++)
    *f += st::real(raw_tmp[i*lda_tmp+i]);
  *f = mt::sqrt(*f);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(tmp,iflag),*iflag);
}

