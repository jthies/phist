/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
static void SUBR(orthogrrfused_cholrr)(TYPE(sdMat_ptr) RR, TYPE(sdMat_ptr) R_1, TYPE(sdMat_ptr) WR_1, 
        int* rank, _MT_ rankTol, int* iflag)
{
#include "phist_std_typedefs.hpp"
  int m = 0, n = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(RR,&m,iflag),*iflag);
  // create R and R_1 from RR
  *rank = 0;
  {
    // stable rank-revealing cholesky
    int perm[m];
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(RR,perm,rank,rankTol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(RR,perm,*rank,R_1,iflag),*iflag);

    if( WR_1 != NULL )
    {
      TYPE(sdMat_ptr) WR_1t = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(WR_1,&n,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&WR_1t,m,n,NULL,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMatT_add_sdMat)(-st::one(),WR_1,st::zero(),WR_1t,iflag),*iflag);
      //PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(RR,perm,*rank,WR_1t,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_forwardSubst_sdMat)(RR,perm,*rank,WR_1t,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMatT_add_sdMat)(st::one(),WR_1t,st::zero(),WR_1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(WR_1t,iflag),*iflag);
    }
  }
//PHIST_CHK_IERR(SUBR(sdMat_print)(R_1,iflag),*iflag);
}


// Q*R_1 = V - W*R_2 with R_2=W'*V
// result returned in place as V
// correct V'V must be supplied and is returned!
//
// Trying to exploit one-step kernel:
// V <- W*M + V*N
// from
// (V-W*WtV)'*(V-W*WtV) = VtV - VtW*WtV - VtW*WtV + VtW*WtW*WtV = VtV - 2 VtW*WtV + VtW*WtV = VtV - VtW*WtV
// for previously orthog. V (wrt. itself) -> I - WtV'*WtV
void SUBR(orthogrrfused)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, 
        TYPE(const_sdMat_ptr) VtV, _MT_ rankTol, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    // get dimensions
    int m, k;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
    phist_const_comm_ptr comm;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);
    // create matrices
    TYPE(sdMat_ptr) R = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R,m,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) R_1 = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R_1,m,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) WtV = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&WtV,k,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) WtV_ = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&WtV_,k,m,comm,iflag),*iflag);

    // calculate first R factor
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
    int Vrank = 0;
    PHIST_CHK_IERR(SUBR(orthogrrfused_cholrr)(R,R_1,NULL,&Vrank,rankTol,iflag),*iflag);

    // calculate WtV, orthog V
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(fused_mvsdi_mvTmv)(st::one(),W,V,R_1,st::zero(),WtV,iflag),*iflag);

    // calculate R2 = W'*V = WtV*R
    PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtV,R,st::zero(),R2,iflag),*iflag);
    // Calculate R1 = R
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R,st::zero(),R1,iflag),*iflag);

    // calculate R factor
    // R <- I - WtV'WtV
    int WVrank = 0;
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(-st::one(),WtV,WtV,st::one(),R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),WtV,st::zero(),WtV_,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(orthogrrfused_cholrr)(R,R_1,WtV_,&WVrank,rankTol,iflag),*iflag);

    // respect R_1 and sign in WtV
    //PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(-st::one(),WtV,R_1,st::zero(),WtV_,iflag),*iflag);

    // all in one fused kernel
    // V <- V*R_1 + W*WtV_
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(W,WtV_,V,R_1,iflag),*iflag);


    // update R1 = R * R1
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R_1,st::zero(),R1,iflag),*iflag);


    PHIST_CHK_IERR(SUBR(sdMat_delete)(R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(WtV,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(WtV_,iflag),*iflag);

    // the return value of this function is the rank of the null space of [W V] on entry
    *iflag=m-WVrank;
}
