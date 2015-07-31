static void SUBR(orthogrr_cholrr)(TYPE(sdMat_ptr) RR, TYPE(sdMat_ptr) R_1, int* rank, int* iflag)
{
  int m = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(RR,&m,iflag),*iflag);
  // create R and R_1 from RR
  *rank = 0;
  {
    // stable rank-revealing cholesky
    int perm[m];
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(RR,perm,rank,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(RR,perm,*rank,R_1,iflag),*iflag);
  }
PHIST_CHK_IERR(SUBR(sdMat_print)(R_1,iflag),*iflag);
}


// Q*R_1 = V - W*R_2 with R_2=W'*V
// result returned in place as V
// correct V'V must be supplied and is returned!
void SUBR(orthogrr)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, TYPE(sdMat_ptr) VtV, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    // get dimensions
    int m, k;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
    const_comm_ptr_t comm;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);
    // create matrices
    TYPE(sdMat_ptr) R = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R,m,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) R_1 = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R_1,m,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) WtV = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&WtV,k,m,comm,iflag),*iflag);

    // we already have the current VtV!
    // so calculate first "R" factor
    int Vrank = 0;
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&Vrank,iflag),*iflag);

    // V <- V*R_1 and WtV <- W'*V
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace)(st::one(),W,V,R_1,st::zero(),WtV,iflag),*iflag);

    // calculate R2 = W'*V = WtV*R
    PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtV,R,st::zero(),R2,iflag),*iflag);
    // calculate R1 = (I- WtV'*WtV)*R = R - WtV'*R2
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R,st::zero(),R1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(-st::one(),WtV,R2,st::one(),R1,iflag),*iflag);

    // V <- V - W*WtV, updating VtV
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_augmented)(-st::one(),W,WtV,st::one(),V,VtV,iflag),*iflag);

    // calculate new R factor
    int WVrank = 0;
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&WVrank,iflag),*iflag);

    // TODO: possibly perform another full sweep with W, if R_1 is too badly conditioned

    // TODO: possibly omit "scaling" step, if VtV is already ok
    // V <- V*R_1, VtV <- V'*V
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace)(st::one(),V,V,R_1,st::zero(),VtV,iflag),*iflag);
    // update R1
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(),R_1,R,st::zero(),R1,iflag),*iflag);


    PHIST_CHK_IERR(SUBR(sdMat_delete)(R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(WtV,iflag),*iflag);

    // the return value of this function is the rank of the null space of V on entry
    *iflag=m-WVrank;
}
