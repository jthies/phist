/*! \file prec_kernels_def.hpp
 * Generic Implementation of interface in phist_prec_kernels.h
 * (should eventually be used by e.g. GHOST as well).
 * 
 * \author "Jonas Thies" <Jonas.Thies@DLR.de>
 * This implementation can be used by any kernel lib that provides
 * PHIST_HIGH_PRECISION_KERNELS and implements sdMat_extract_err.
*/


////////////////////////////////////////////////////////////////////////////////////////////////
// implementation of public interface to kernels in prec_kernels.c                            //
////////////////////////////////////////////////////////////////////////////////////////////////

void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) M, int* perm, int* rank, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  lidx_t ldM, n,m;
  _ST_ *Mval, *Merr;
  *iflag=0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(M,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&Mval,&ldM,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error)(M,&Merr,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(prec_cholesky)(Mval,Merr,m,ldM,perm,rank,iflag),*iflag);
}

void SUBR(sdMat_backwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  lidx_t ldR, n, m, ldX, k;
  _ST_ *Rval, *Rerr, *Xval, *Xerr;
  *iflag=0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&Rval,&ldR,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error)(R,&Rerr,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(X,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(X,&k,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(X,&Xval,&ldX,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error)(X,&Xerr,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(prec_backwardSubst)(Rval,Rerr,n,ldR,perm,rank,Xval,Xerr,k,ldX,iflag),*iflag);
}

//! forward substitution. \ingroup prec

//! forward substitution for pivoted conj. transposed upper triangular cholesky factor
void SUBR(sdMat_forwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  lidx_t ldR, n, m, ldX, k;
  _ST_ *Rval, *Rerr, *Xval, *Xerr;
  *iflag=0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&Rval,&ldR,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error)(R,&Rerr,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(X,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(X,&k,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(X,&Xval,&ldX,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error)(X,&Xerr,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(prec_forwardSubst)(Rval,Rerr,n,ldR,perm,rank,Xval,Xerr,k,ldX,iflag),*iflag);
}

//! given B=V'V, compute (in place) B^ s.t. V*B^ is orthonormal. The rank of V is returned in *rank.
void SUBR(sdMat_qb)(TYPE(sdMat_ptr) B, 
                    TYPE(sdMat_ptr) B_1, 
                    int* rank, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  lidx_t ldB, ldB_1, n, m;
  _ST_ *Bval, *B_1val, *Berr, *B_1err;
  *iflag=0;
  
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(B,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(B,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B,&Bval,&ldB,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error)(B,&Berr,iflag),*iflag);
  if (B_1!=NULL)
  {
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B_1,&B_1val,&ldB_1,iflag),*iflag);
    PHIST_CHK_IERR(*iflag=(ldB==ldB_1)?0:-1,*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_error)(B_1,&B_1err,iflag),*iflag);
  }
  else
  {
    // B^{-1} not wanted
    B_1val=NULL;
    B_1err=NULL;
  }
  PHIST_CHK_IERR(SUBR(prec_qb)(Bval,Berr,B_1val,B_1err,n,ldB,rank,iflag),*iflag);
}
