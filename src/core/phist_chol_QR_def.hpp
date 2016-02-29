// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
// with rank-revealing pivoted cholesky (SVRR)
// If the input vector has m columns and rank r, iflag=m-r is
// returned. The nullspace of V is randomized and orthogonalized
// against the other columns. If the output matrix is denoted
// by Q, Q and V are related by Q*R=V.
//
// This routine is the fallback kernel used by orthog if the kernel
// library doesn't provide mvec_QR.
void SUBR(chol_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    int m;
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    const_comm_ptr_t comm;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);
    TYPE(sdMat_ptr) R_1 = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R_1,m,m,comm,iflag),*iflag);

    // S=V'V
    if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),R,iflag),*iflag);
    int perm[m];
    int rank;
    if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(R,perm,&rank,iflag),*iflag);

    // construct inv(R)
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(R,perm,rank,R_1,iflag),*iflag);
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,R_1,iflag),*iflag);
    int newRank = rank;
    int nIter = 0;
    while(newRank < m && nIter++ < 2)
    {
      TYPE(mvec_ptr) V_=NULL;
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&V_,newRank,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_random)(V_,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_delete)(V_,iflag),*iflag);
      TYPE(sdMat_ptr) R_=NULL;
      PHIST_CHK_IERR(SUBR(sdMat_create)(&R_,m,m,comm,iflag),*iflag);
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),R_,iflag),*iflag);
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(sdMat_cholesky)(R_,perm,&newRank,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(R_,perm,newRank,R_1,iflag),*iflag);
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,R_1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R_,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    *iflag=m-rank;
    return;
}
