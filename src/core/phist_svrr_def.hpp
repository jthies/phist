// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
// with rank-revealing pivoted cholesky (SVRR)
// If the input vector has m columns and rank r, iflag=m-r is
// returned. Columns 0:iflag-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. The third argument, E, 
// should be preallocated by the user with m elements, m being
// the number of columns in V. On successful return (iflag>=0),
// e[j] indicates the norm of V(:,j) before the orthogonali-  
// zation step. On exit, V(:,j), j>*iflag has 2-norm 1.
void SUBR(svrr)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    int m;
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    phist_const_comm_ptr comm;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);
    TYPE(sdMat_ptr) R = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R,m,m,comm,iflag),*iflag);

    // S=V'V
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),R,iflag),*iflag);
#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_SOUT(PHIST_INFO,"Q^TQ:\n");
PHIST_CHK_IERR(SUBR(sdMat_print)(R,iflag),*iflag);
PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R,st::zero(),B,iflag),*iflag);
#endif

    // create B from R
    int rank = 0;
    {
      // stable rank-revealing cholesky
      int perm[m];
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(sdMat_cholesky)(R,perm,&rank,iflag),*iflag);
#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_SOUT(PHIST_INFO,"R^T:\n");
PHIST_CHK_IERR(SUBR(sdMat_print)(R,iflag),*iflag);
PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(st::one(),R,R,-st::one(),B,iflag),*iflag);
PHIST_SOUT(PHIST_INFO,"Q^TQ-RR^T:\n");
PHIST_CHK_IERR(SUBR(sdMat_print)(B,iflag),*iflag);
#endif
      PHIST_CHK_IERR(SUBR(sdMat_identity)(B,iflag),*iflag);
      if( robust )
        *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(R,perm,rank,B,iflag),*iflag);
#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_SOUT(PHIST_INFO,"R^-T:\n");
PHIST_CHK_IERR(SUBR(sdMat_print)(B,iflag),*iflag);
#endif
    }

    PHIST_CHK_IERR(SUBR(sdMat_delete)(R,iflag),*iflag);

    // compute V <- V*B to get an orthogonal V (up to the first (m-rank) columns,
    // which will be exactly zero)
    if( robust )
      *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,B,iflag),*iflag);

// the return value of this function is the rank of the null space of V on entry
*iflag=m-rank;
}
