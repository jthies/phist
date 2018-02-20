/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182)
// with rank-revealing pivoted cholesky (SVRR)
// If the input vector has m columns and rank r, iflag=m-r is
// returned. The nullspace of V is randomized and orthogonalized
// against the other columns. If the output matrix is denoted
// by Q, Q and V are related by Q*R=V. R(:,perm) is upper triangular.
//
// This routine is the fallback kernel used by orthog if the kernel
// library doesn't provide mvec_QR.
extern "C" void SUBR(chol_QRp)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int perm[], int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    int m;
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    MT rankTol = mt::rankTol(robust);
    *iflag=0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    phist_const_comm_ptr comm;
    
    for (int i=0; i<m; i++) perm[i]=i;
    
    if (m==1)
    {
      // just normalize the vector
      MT nrm;
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvec_normalize)(V,&nrm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,(ST)nrm,iflag),*iflag);
      int rank=1;
      if (nrm<rankTol)
      {
        // randomize the vector
        PHIST_CHK_IERR(SUBR(mvec_random)(V,iflag),*iflag);
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
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

    // S=V'V
    if (robust) *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),R,iflag),*iflag);

    int rank;
    if (robust) *iflag = PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(R,perm,&rank,rankTol,iflag),*iflag);
    // since R is output, sync it with devices
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(R,iflag),*iflag);

    // construct inv(R)
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(R,perm,rank,R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(R_1,iflag),*iflag);

    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,R_1,iflag),*iflag);

    if (rank < m)
    {
      // randomize the null space
      TYPE(mvec_ptr) Vorth=NULL,V_=NULL;
      if (rank)
      {
        PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&Vorth,0,rank-1,iflag),*iflag);
      }
      PHIST_CHK_IERR(SUBR(mvec_view_block)(V,&V_,rank,m-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_random)(V_,iflag),*iflag);
      // orthogonalize it by recursively calling this function
      TYPE(sdMat_ptr) R1=NULL, R2=NULL;
      if (rank)
      {
        PHIST_CHK_IERR(SUBR(sdMat_create)(&R2,rank,m-rank,comm,iflag),*iflag);
      }
      PHIST_CHK_IERR(SUBR(sdMat_create)(&R1,m-rank,m-rank,comm,iflag),*iflag);

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
        SUBR(chol_QR)(V_,R1,iflag);
        if (*iflag)
        {
          PHIST_SOUT(PHIST_ERROR,"catastrophic breakdown in chol_QR, returning -1\n");
          *iflag=-1;
          return;
        }
        if (rank>0)
        {
          if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
          PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),Vorth,V_,st::zero(),R2,iflag),*iflag);
          if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
          PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Vorth,R2,st::one(),V_,iflag),*iflag);
        }
      }
      // The orthogonalization coefficients are simply thrown away because these columns 
      // do not contribute to recovering V by V=Q*R
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R1,iflag),*iflag);
      if (rank>0)
      {
        PHIST_CHK_IERR(SUBR(sdMat_delete)(R2,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_delete)(Vorth,iflag),*iflag);
      }
      PHIST_CHK_IERR(SUBR(mvec_delete)(V_,iflag),*iflag);
    }
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    *iflag=m-rank;
    return;
}

//! this variant of chol_QR discards the perm array. It has the same interface
//! as mvec_QR but returns a column-permuted upper triangular matrix R.
void SUBR(chol_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  int m;
  int iflag_in=*iflag;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
  int perm[m];
  *iflag=iflag_in;
  SUBR(chol_QRp)(V,R,perm,iflag);
}
