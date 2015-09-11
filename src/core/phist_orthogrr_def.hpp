
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
void SUBR(orthogrr)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, TYPE(const_sdMat_ptr) WtW_I, TYPE(sdMat_ptr) VtV, _MT_ desiredEps, int maxIter, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    // get dimensions
    int m, k = 0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    if( W != NULL ) {PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);}
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
    if( k > 0 ) {PHIST_CHK_IERR(SUBR(sdMat_create)(&WtV,k,m,comm,iflag),*iflag);}
    TYPE(sdMat_ptr) VtV_I = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&VtV_I,m,m,comm,iflag),*iflag);
    // if we have a highly accurate approximation of the orthogonalization error of W
    // (e.g. WtW_I), we can project it out to obtain a smaller error WtV
    TYPE(sdMat_ptr) WtW_inv = NULL;
    TYPE(sdMat_ptr) EWtV = NULL;
    int rankWtW, *permWtW = NULL;
    if( k > 0 && WtW_I )
    {
      PHIST_CHK_IERR(SUBR(sdMat_create)(&WtW_inv,k,k,comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&EWtV,k,m,comm,iflag),*iflag);
      // construct cholesky factorization of WtW
      PHIST_CHK_IERR(SUBR(sdMat_identity)(WtW_inv,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),WtW_I,st::one(),WtW_inv,iflag),*iflag);
      permWtW = new int[k];
      PHIST_CHK_IERR(SUBR(sdMat_cholesky)(WtW_inv,permWtW,&rankWtW,iflag),*iflag);
    }

    int rank = m;
    MT VtV_err = mt::zero();
    MT WtV_err = mt::zero();
    PHIST_CHK_IERR(SUBR(sdMat_identity)(VtV_I,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,-st::one(),VtV_I,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_normF)(VtV_I,&VtV_err,iflag),*iflag);

    // initialize R1 and R2
    if( R1 != NULL ) {PHIST_CHK_IERR(SUBR(sdMat_identity)(R1,iflag),*iflag);}
    if( R2 != NULL ) {PHIST_CHK_IERR(SUBR(sdMat_put_value)(R2,st::zero(),iflag),*iflag);}

    // allow multiple sweeps (2 should be enough if high prec is used!)
    int iter = 0;
    for(; iter < maxIter; iter++)
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
      if( k == 0 )
      {
        iter++;
        break;
      }

      // V <- V*R_1 and WtV <- W'*V
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace)(st::one(),W,V,R_1,st::zero(),WtV,iflag),*iflag);
      if( WtW_I )
      {
        // correct WtV by error of WtW from
        // 0 = W^T(V-W(WtV+X))
        //   = WtV - WtW (WtV+X)
        //   = WtV - (I+E) (WtV+X)
        //   = -E WtV - (I+E)X
        // So we obtain the correction
        // X = - WtW^(-1) WtW_I WtV
 
        PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtW_I,WtV,st::zero(),EWtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(WtW_inv,permWtW,rankWtW,EWtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_forwardSubst_sdMat)(WtW_inv,permWtW,rankWtW,EWtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(-st::one(),EWtV,st::one(),WtV,iflag),*iflag);
      }

      // update R1 <- R * R1
      if( R1 != NULL )
      {
        PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1_tmp,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R1_tmp,st::zero(),R1,iflag),*iflag);
      }
      // check if we really need another orthogonalization step with W
      VtV_err = mt::eps();
      PHIST_CHK_IERR(SUBR(sdMat_normF)(WtV,&WtV_err,iflag),*iflag);
      PHIST_SOUT(PHIST_INFO, "orthogRR: iter %d phase 1, desired eps %8.4e, WtV err. %8.4e, (est.) VtV err. %8.4e\n", iter, desiredEps, WtV_err, VtV_err);
      if( WtV_err <= desiredEps )
      {
        iter++;
        break;
      }

      // update R2 <- R2 + W'*V = R2 + WtV*R
      if( R2 != NULL ) {PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtV,R,st::one(),R2,iflag),*iflag);}

      // V <- V - W*WtV, updating VtV
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat_augmented)(-st::one(),W,WtV,st::one(),V,VtV,iflag),*iflag);

      // calculate new R factor
      int WVrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
#ifdef ORTHOGRR_USE_SVQB
      PHIST_CHK_IERR(SUBR(orthogrr_svqb)(R,R_1,&WVrank,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&WVrank,iflag),*iflag);
#endif
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
      {
        iter++;
        break;
      }
    }

    // continue iteration with VtV
    for(; iter <= maxIter; iter++)
    {
      if( VtV_err <= desiredEps )
        break;

      // V <- V*R_1, VtV <- V'*V
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace)(st::one(),V,V,R_1,st::zero(),VtV,iflag),*iflag);

      // update R1 <- R * R1
      if( R1 != NULL )
      {
        PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1_tmp,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R1_tmp,st::zero(),R1,iflag),*iflag);
      }

      // calculate new R factor
      int Vrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
#ifdef ORTHOGRR_USE_SVQB
      PHIST_CHK_IERR(SUBR(orthogrr_svqb)(R,R_1,&Vrank,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&Vrank,iflag),*iflag);
#endif
      rank = std::min(rank,Vrank);

      // calculate error
      PHIST_CHK_IERR(SUBR(sdMat_identity)(VtV_I,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,-st::one(),VtV_I,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_normF)(VtV_I,&VtV_err,iflag),*iflag);
      PHIST_SOUT(PHIST_INFO, "orthogRR: iter %d phase 3, desired eps %8.4e, VtV err. %8.4e\n", iter-1, desiredEps, VtV_err);
    }


    PHIST_CHK_IERR(SUBR(sdMat_delete)(R,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_tmp,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    if( WtV )     {PHIST_CHK_IERR(SUBR(sdMat_delete)(WtV,    iflag),*iflag);}
    if( EWtV )    {PHIST_CHK_IERR(SUBR(sdMat_delete)(EWtV,   iflag),*iflag);}
    if( WtW_inv ) {PHIST_CHK_IERR(SUBR(sdMat_delete)(WtW_inv,iflag),*iflag);}
    PHIST_CHK_IERR(SUBR(sdMat_delete)(VtV_I,iflag),*iflag);

    // the return value of this function is the rank of the null space of V on entry
    *iflag=m-rank;
}
