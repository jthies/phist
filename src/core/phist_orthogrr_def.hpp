
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
//PHIST_CHK_IERR(SUBR(sdMat_print)(R_1,iflag),*iflag);
}

static void SUBR(orthogrr_svqb)(TYPE(sdMat_ptr) RR, TYPE(sdMat_ptr) R_1, int* rank, int* iflag)
{
#include "phist_std_typedefs.hpp"
  // compute B s.t. V*B=Q is orthogonal, with RR=V'V on input
  // and inv(B) on output, and R_1=B on output.
  
  // we first copy the input matrix because B and B_1 are exchanged in the definition of
  // the kernel routine sdMat_qb:
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),RR,st::zero(),R_1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_qb)(R_1,RR,rank,iflag),*iflag);
//PHIST_CHK_IERR(SUBR(sdMat_print)(R_1,iflag),*iflag);
}


// Q*R_1 = V - W*R_2 with R_2=W'*V
// result returned in place as V
// correct V'V must be supplied and is returned!
void SUBR(orthogrr)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, TYPE(sdMat_ptr) VtV, _MT_ desiredEps, int maxIter, int* iflag)
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
    TYPE(sdMat_ptr) R1_tmp = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1_tmp,m,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) R_1 = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R_1,m,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) WtV = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&WtV,k,m,comm,iflag),*iflag);
    TYPE(sdMat_ptr) VtV_I = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&VtV_I,m,m,comm,iflag),*iflag);

    int rank = m;
    MT VtV_err = mt::zero();
    MT WtV_err = mt::zero();

    // initialize R1 and R2
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R2,st::zero(),iflag),*iflag);

    // allow multiple sweeps (2 should be enough if high prec is used!)
    for(int iter = 0; iter < maxIter; iter++)
    {
      // we already have the current VtV!
      // so calculate first "R" factor
      int Vrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
#ifdef ORTHOGRR_USE_SVQB
      PHIST_CHK_IERR(SUBR(orthogrr_svqb)(R,R_1,&Vrank,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&Vrank,iflag),*iflag);
#endif
      rank = std::min(rank,Vrank);

      // V <- V*R_1 and WtV <- W'*V
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace)(st::one(),W,V,R_1,st::zero(),WtV,iflag),*iflag);

      // update R1 <- R * R1
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1_tmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R1_tmp,st::zero(),R1,iflag),*iflag);
      // check if we really need another orthogonalization step with W
      VtV_err = mt::eps();
      PHIST_CHK_IERR(SUBR(sdMat_normF)(WtV,&WtV_err,iflag),*iflag);
      PHIST_SOUT(PHIST_INFO, "orthogRR: iter %d phase 1, desired eps %8.4e, WtV err. %8.4e, (est.) VtV err. %8.4e\n", iter, desiredEps, WtV_err, VtV_err);
      if( WtV_err <= desiredEps )
        break;

      // update R2 <- R2 + W'*V = R2 + WtV*R
      PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtV,R,st::one(),R2,iflag),*iflag);

      // V <- V - W*WtV, updating VtV
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat_augmented)(-st::one(),W,WtV,st::one(),V,VtV,iflag),*iflag);

      // calculate new R factor
      int WVrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&WVrank,iflag),*iflag);
      rank = std::min(rank,WVrank);

      // estimate error
      PHIST_CHK_IERR(SUBR(sdMat_normF)(R_1,&WtV_err,iflag),*iflag);
      WtV_err *= mt::eps();
      PHIST_CHK_IERR(SUBR(sdMat_identity)(VtV_I,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,-st::one(),VtV_I,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_normF)(VtV_I,&VtV_err,iflag),*iflag);

      // possibly perform another full sweep with W, if R_1 is too badly conditioned
      PHIST_SOUT(PHIST_INFO, "orthogRR: iter %d phase 2, desired eps %8.4e, (est.) WtV err. %8.4e, VtV err. %8.4e\n", iter, desiredEps, WtV_err, VtV_err);
      if( WtV_err <= desiredEps )
        break;
    }

    // possibly omit "scaling" step, if VtV is already ok
    if( VtV_err > desiredEps )
    {
      // V <- V*R_1, VtV <- V'*V
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace)(st::one(),V,V,R_1,st::zero(),VtV,iflag),*iflag);

      // update R1 <- R * R1
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1_tmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R1_tmp,st::zero(),R1,iflag),*iflag);
    }


    PHIST_CHK_IERR(SUBR(sdMat_delete)(R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_tmp,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(WtV,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(VtV_I,iflag),*iflag);

    // the return value of this function is the rank of the null space of V on entry
    *iflag=m-rank;
}
