// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
// with rank-revealing pivoted cholesky (SVRR)
// If the input vector has m columns and rank r, iflag=m-r is
// returned. The nullspace of V is randomized and orthogonalized
// against the other columns. If the output matrix is denoted
// by Q, Q and V are related by Q*R=V.
//
// This routine is the fallback kernel used by orthog if the kernel
// library doesn't provide mvec_QR.
extern "C" void SUBR(chol_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    int m;
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    const_comm_ptr_t comm;
    
    if (m==1)
    {
      // just normalize the vector
      MT nrm;
      PHIST_CHK_IERR(SUBR(mvec_normalize)(V,&nrm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,(ST)nrm,iflag),*iflag);
      int rank=1;
      if (nrm<10*mt::eps())
      {
        // randomize the vector
        PHIST_CHK_IERR(SUBR(mvec_random)(V,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_normalize)(V,&nrm,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),iflag),*iflag);
        rank=0;
      }
      *iflag=1-rank;
      return;
    }
    
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);
    TYPE(sdMat_ptr) R_1 = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R_1,m,m,comm,iflag),*iflag);
#if PHIST_OUTLEV>=PHIST_DEBUG
int nrows;
PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nrows,iflag),*iflag);
if (nrows<=100)
{
PHIST_SOUT(PHIST_INFO,"%%CHOLQR\nV=[...\n");
SUBR(mvec_print)(V,iflag);
PHIST_SOUT(PHIST_INFO,"];\n");
}
#endif
    // S=V'V
    if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),R,iflag),*iflag);
    int perm[m];
    int rank;
    if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(R,perm,&rank,iflag),*iflag);
#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_SOUT(PHIST_INFO,"%%CHOLQR: cholesky(V'V) has rank %d\nR=[",rank);
SUBR(sdMat_print)(R,iflag);
PHIST_SOUT(PHIST_INFO,"];\n");
#endif
    // construct inv(R)
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(R,perm,rank,R_1,iflag),*iflag);
#if PHIST_OUTLEV>=PHIST_DEBUG
PHIST_SOUT(PHIST_INFO,"%%CHOLQR: inv(R)=\nR_1=[");
SUBR(sdMat_print)(R_1,iflag);
PHIST_SOUT(PHIST_INFO,"];\n");
#endif
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,R_1,iflag),*iflag);

#if PHIST_OUTLEV>=PHIST_DEBUG
if (nrows<=100)
{
PHIST_SOUT(PHIST_INFO,"CHOLQR: Q=V*R_I (before randomizing null space)\nQ0=[...");
SUBR(mvec_print)(V,iflag);
PHIST_SOUT(PHIST_INFO,"];\n");
}
#endif
    int newRank = rank;
    int nIter = 0;
    while(newRank < m && nIter++ < 2)
    {
      // randomize the null space
      TYPE(mvec_ptr) Vorth=NULL,V_=NULL;
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vorth,0,newRank-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&V_,newRank,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_random)(V_,iflag),*iflag);
      // orthogonalize it by recursively calling this function
      TYPE(sdMat_ptr) R1=NULL, R2=NULL;
      PHIST_CHK_IERR(SUBR(sdMat_create)(&R1,m-newRank,m-newRank,comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&R2,newRank,m-newRank,comm,iflag),*iflag);

      // normalize the new random vectors and project out the existing orthogonal columns.
      // We just do two classical Gram-Schmidt+CholQR steps here, assuming that this is a fallback kernel
      // if the kernel lib doesn't provide TSQR or high precision and fused kernels (TODO: this can be 
      // done in smarter ways, see orthogrr and orthogrrfused). 
      for (int ii=0;ii<2;ii++)
      {
        // This call will report an error if the random columns created are
        // not linearly independent. We assume that this will not happen in
        // practice
        if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(SUBR(chol_QR)(V_,R1,iflag),*iflag);
        if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vorth,V_,st::zero(),R2,iflag),*iflag);
        if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Vorth,R2,st::one(),V_,iflag),*iflag);
      }
      // The orthogonalization coefficients are simply thrown away because these columns 
      // do not contribute to recovering V by V=Q*R
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R2,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    *iflag=m-rank;
    return;
}
