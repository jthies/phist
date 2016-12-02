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

extern "C" void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) M, int* perm, int* rank, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  phist_lidx ldM, n,m;
  _ST_ *Mval, *Merr=NULL;
#if PHIST_HIGH_PRECISION_KERNELS_FORCE
  bool robust=true;
#else
  bool robust=(*iflag&PHIST_ROBUST_REDUCTIONS);
#endif
  *iflag=0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(M,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:PHIST_INVALID_INPUT,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&Mval,&ldM,iflag),*iflag);
#ifdef PHIST_HIGH_PRECISION_KERNELS
  SUBR(sdMat_extract_error)(M,&Merr,iflag);
  if (robust&&Merr!=NULL)
  {
    PHIST_CHK_IERR(SUBR(prec_cholesky)(Mval,Merr,m,ldM,perm,rank,iflag),*iflag);
  }
  else
#endif
  {
    PHIST_CHK_IERR(SUBR(cholesky)(Mval,m,ldM,perm,rank,iflag),*iflag);
  }
}

extern "C" void SUBR(sdMat_backwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  phist_lidx ldR, n, m, ldX, k;
  _ST_ *Rval, *Rerr=NULL, *Xval, *Xerr=NULL;
#if PHIST_HIGH_PRECISION_KERNELS_FORCE
  bool robust=true;
#else
  bool robust=(*iflag&PHIST_ROBUST_REDUCTIONS);
#endif
  *iflag=0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&Rval,&ldR,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(X,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(X,&k,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(X,&Xval,&ldX,iflag),*iflag);
#ifdef PHIST_HIGH_PRECISION_KERNELS
  SUBR(sdMat_extract_error)(R,&Rerr,iflag);
  SUBR(sdMat_extract_error)(X,&Xerr,iflag);

  if (robust&&(Rerr!=NULL)&&(Xerr!=NULL))
  {
    PHIST_CHK_IERR(SUBR(prec_backwardSubst)(Rval,Rerr,n,ldR,perm,rank,Xval,Xerr,k,ldX,iflag),*iflag);
  }
  else
#endif
  {
    PHIST_CHK_IERR(SUBR(backwardSubst)(Rval,n,ldR,perm,rank,Xval,k,ldX,iflag),*iflag);
  }
}

//! forward substitution. \ingroup prec

//! forward substitution for pivoted conj. transposed upper triangular cholesky factor
extern "C" void SUBR(sdMat_forwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  phist_lidx ldR, n, m, ldX, k;
  _ST_ *Rval, *Rerr=NULL, *Xval, *Xerr=NULL;
#if PHIST_HIGH_PRECISION_KERNELS_FORCE
  bool robust=true;
#else
  bool robust=(*iflag&PHIST_ROBUST_REDUCTIONS);
#endif
  *iflag=0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&Rval,&ldR,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(X,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(X,&k,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(X,&Xval,&ldX,iflag),*iflag);

#ifdef PHIST_HIGH_PRECISION_KERNELS
  SUBR(sdMat_extract_error)(R,&Rerr,iflag);
  SUBR(sdMat_extract_error)(X,&Xerr,iflag);

  if (robust&&(Rerr!=NULL)&&(Xerr!=NULL))
  {
    PHIST_CHK_IERR(SUBR(prec_forwardSubst)(Rval,Rerr,n,ldR,perm,rank,Xval,Xerr,k,ldX,iflag),*iflag);
  }
  else
#endif
  {
    PHIST_CHK_IERR(SUBR(forwardSubst)(Rval,n,ldR,perm,rank,Xval,k,ldX,iflag),*iflag);
  }
}

//! given B=V'V, compute (in place) B^ s.t. V*B^ is orthonormal. The rank of V is returned in *rank.
extern "C" void SUBR(sdMat_qb)(TYPE(sdMat_ptr) B, 
                    TYPE(sdMat_ptr) B_1, 
                    int* rank, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  phist_lidx ldB, ldB_1, n, m;
  _ST_ *Bval, *B_1val, *Berr=NULL, *B_1err=NULL;
  PHIST_TOUCH(Berr);
  PHIST_TOUCH(B_1err);
#if PHIST_HIGH_PRECISION_KERNELS_FORCE
  bool robust=true;
#else
  bool robust=(*iflag&PHIST_ROBUST_REDUCTIONS);
#endif
  *iflag=0;
  
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(B,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(B,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=(n==m)?0:-1,*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B,&Bval,&ldB,iflag),*iflag);
#ifdef PHIST_HIGH_PRECISION_KERNELS
  SUBR(sdMat_extract_error)(B,&Berr,iflag);
#endif
  if (B_1!=NULL)
  {
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(B_1,&B_1val,&ldB_1,iflag),*iflag);
    PHIST_CHK_IERR(*iflag=(ldB==ldB_1)?0:-1,*iflag);
#ifdef PHIST_HIGH_PRECISION_KERNELS
    SUBR(sdMat_extract_error)(B_1,&B_1err,iflag);
#endif
  }
  else
  {
    // B^{-1} not wanted
    B_1val=NULL;
    B_1err=NULL;
  }
#ifdef PHIST_HIGH_PRECISION_KERNELS
  if (robust&&(Berr!=NULL))
  {
    PHIST_CHK_IERR(SUBR(prec_qb)(Bval,Berr,B_1val,B_1err,n,ldB,rank,iflag),*iflag);
  }
  else
#endif
  {
    PHIST_CHK_IERR(SUBR(qb)(Bval,B_1val,n,ldB,rank,iflag),*iflag);
  }
}

//! computes in-place the inverse of a Hermitian and positive semi-definite matrix using Cholesky factorization.
//! If A is singular (actually semi-definite, that is), the pseudo-inverse is computed using rank-revealing Cholesky.
//! The rank of A on input is returned as *rank.
void SUBR(sdMat_inverse)(TYPE(sdMat_ptr) A_hpd, int* rank, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  int iflag_in=*iflag;
  int n,m;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(A_hpd,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(A_hpd,&m,iflag),*iflag);
  PHIST_CHK_IERR(*iflag = n==m? 0: PHIST_INVALID_INPUT, *iflag);

  // create matrix with NULL communicator
  phist_const_comm_ptr comm=NULL;
  
  TYPE(sdMat_ptr) RR=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&RR,m,m,comm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_set_block)(RR,A_hpd,0,m-1,0,m-1,iflag),*iflag);
  
  int perm[m];
  *iflag=iflag_in;
  PHIST_CHK_IERR(SUBR(sdMat_cholesky)(RR,perm,rank,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_identity)(A_hpd,iflag),*iflag);
  *iflag=iflag_in;
  // compute A^{-1} = (R'R)\I= R'\ R\ I
  PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(RR,perm,*rank,A_hpd,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_forwardSubst_sdMat)(RR,perm,*rank,A_hpd,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(RR,iflag),*iflag);
}

//! computes in-place the (transpose of the Moore-Penrose)
//! pseudo-inverse of an arbitrary m x n matrix. The rank
//! of A on input is returned as *rank.
void SUBR(sdMat_pseudo_inverse)(TYPE(sdMat_ptr) A_gen, int* rank, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  int iflag_in=*iflag;
  bool high_prec = *iflag&PHIST_ROBUST_REDUCTIONS;
  int n,m;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(A_gen,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(A_gen,&n,iflag),*iflag);
  
  phist_const_comm_ptr comm=NULL;
  
  *rank=-1;
  
  TYPE(sdMat_ptr) U=NULL, Sigma=NULL, Vt=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&U,m,m,comm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&Sigma,m,n,comm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&Vt,n,n,comm,iflag),*iflag);
  SdMatOwner<_ST_> _Sigma(Sigma),_U(U),_Vt(Vt);

  // eigenvalue decomposition, A = V*Sigma*W'
  *iflag=iflag_in;
  PHIST_CHK_IERR(SUBR(sdMat_svd)(A_gen,U,Sigma,Vt,iflag),*iflag);

  // make tiny singular values exactly 0, invert the others
  _ST_ *Sigma_raw=NULL, *Sigma_err=NULL;
  phist_lidx ldS;

  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(Sigma,&Sigma_raw,&ldS,iflag),*iflag);
  _MT_ sval_max = st::abs(Sigma_raw[0]), sval_max_err=mt::zero();
  *rank=std::min(m,n);
  if (high_prec)
  {
#ifdef PHIST_HIGH_PRECISION_KERNELS
    PHIST_CHK_IERR(SUBR(sdMat_extract_error)(Sigma,&Sigma_err,iflag),*iflag);
    sval_max_err = st::abs(Sigma_err[0]);
#else
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
#endif
  }
  for (int i=0; i<std::min(n,m); i++)
  {
    _ST_ sval = Sigma_raw[i*ldS+i];
#ifdef PHIST_HIGH_PRECISION_KERNELS
    if (high_prec)
    {
      _ST_ serr = Sigma_err[i*ldS+i];
      if (st::abs(sval)<sval_max*mt::rankTol(high_prec))
      {
        Sigma_raw[i*ldS+i]=st::zero();
        Sigma_err[i*ldS+i]=st::zero();
        (*rank)--;
      }
      else
      {
        // compoute 1/sigma_i
        _ST_ s,sC;
        st::prec_div(Sigma_raw[i*ldS+i], Sigma_err[i*ldS+i],s,sC);
        Sigma_raw[i*ldS+i]=s;
        Sigma_err[i*ldS+i]=sC;
      }
    }
    else
#endif
    {
      if (st::abs(sval)<sval_max*mt::rankTol())
      {
        Sigma_raw[i*ldS+i]=st::zero();
        (*rank)--;
      }
      else
      {
        Sigma_raw[i*ldS+i]=st::one()/sval;
      }
    }
  }
  // multiply the three matrices back together to get the Pseudo-Inverse,
  // B=A^{+,T} <- U*inv(Sigma)*V'
  
  // USig <- U*inv(Sigma)
  TYPE(sdMat_ptr) USig=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&USig,m,n,comm,iflag),*iflag);
  SdMatOwner<_ST_> _USig(USig);
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(Sigma,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(U,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),U,Sigma,st::zero(),USig,iflag),*iflag);
  
  // A <- U*inv(Sigma)*V'
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(Vt,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),USig,Vt,st::zero(),A_gen,iflag),*iflag);
  
  return;  
}

//! singular value decomposition, A = U*Sigma*Vt  

//! A, Sigma are  m x n, U m x m, V n x n.          
//! The function just calls the lapack routine XGESVD with JOBU=JOBV='A',
//! so all the left and right singular vectors are computed and A is filled   
//! with garbage on output. The transpose of V is returned, Vt = transpose(V),
//! and the singular values are sorted on the diagonal of Sigma by decreading
//! magnitude.
//!
//! We also allow Sigma to have dimension (min(m,n),1), in that case we return the
//! diagonal entries of Sigma only (the actual singular values).
void SUBR(sdMat_svd)(TYPE(sdMat_ptr) A, TYPE(sdMat_ptr) U, TYPE(sdMat_ptr) Sigma, TYPE(sdMat_ptr) Vt, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  int iflag_in=*iflag;
  bool high_prec = *iflag&PHIST_ROBUST_REDUCTIONS;
  int n,m;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(A,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(A,&n,iflag),*iflag);  
  
  int nrowsSigma, ncolsSigma;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(Sigma,&nrowsSigma,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(Sigma,&ncolsSigma,iflag),*iflag);  
  bool svals_only = nrowsSigma<m || ncolsSigma<n;
  
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(Sigma,st::zero(),iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(Sigma,iflag),*iflag);
  
  _ST_ *A_val, *U_val, *Vt_val, *S_val;
  int ldA, ldU, ldVt, ldSigma;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(A,&A_val,&ldA,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(Sigma,&S_val,&ldSigma,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(U,&U_val,&ldU,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(Vt,&Vt_val,&ldVt,iflag),*iflag);
  
  if (high_prec)
  {
#ifdef PHIST_HIGH_PRECISION_KERNELS
    _ST_ *A_err, *U_err, *Vt_err, *S_err;
    PHIST_CHK_IERR(SUBR(sdMat_extract_error)(A,&A_err,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_error)(Sigma,&S_err,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_error)(U,&U_err,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_extract_error)(Vt,&Vt_err,iflag),*iflag);
    int mn=std::min(m,n);
    _MT_ RS_val[mn],RS_err[mn];
#if defined(IS_DOUBLE)&&(!defined(IS_COMPLEX))&&defined(PHIST_HAVE_MPACK_QD)
    phist_Drgesvd("A","A",m,n,A_val,A_err,ldA,RS_val,RS_err,U_val,U_err,ldU,Vt_val,Vt_err,ldVt,iflag);
#else
    *iflag=PHIST_NOT_IMPLEMENTED;
#endif
    // copy the returned singular values to the diagonal of Sigma
    int ldS=svals_only?0:ldSigma;
    for (int i=0; i<mn; i++)
    {
      S_val[i*ldS+i]=(_ST_)RS_val[i];
      S_err[i*ldS+i]=(_ST_)RS_err[i];
    }
#else
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
#endif  
  }
  else
  {
    // create work array
    _ST_* work=NULL;
    int lwork=-1;
    _ST_ tmp_work;
    int mn=std::min(m,n);
    _MT_ RS_val[mn];
    const char jobu='A', jobvt='A';
#ifdef IS_COMPLEX
    _MT_ rwork[5*mn];
    PHIST_TG_PREFIX(GESVD)((phist_blas_char*)(&jobu),(phist_blas_char*)(&jobvt),&m,&n,
        (st::blas_scalar_t*)A_val,&ldA,
        (mt::blas_scalar_t*)RS_val,
        (st::blas_scalar_t*)U_val,&ldU,
        (st::blas_scalar_t*)Vt_val,&ldVt,
        (st::blas_scalar_t*)&tmp_work,&lwork,
        (mt::blas_scalar_t*)rwork,iflag);
#else
    PHIST_TG_PREFIX(GESVD)((phist_blas_char*)(&jobu),(phist_blas_char*)(&jobvt),&m,&n,A_val,&ldA,RS_val,U_val,&ldU,Vt_val,&ldVt,&tmp_work,&lwork,iflag);
#endif
    lwork=(int)st::real(tmp_work);
    work=new _ST_[lwork];
#ifdef IS_COMPLEX
    PHIST_TG_PREFIX(GESVD)((phist_blas_char*)(&jobu),(phist_blas_char*)(&jobvt),&m,&n,
        (st::blas_scalar_t*)A_val,&ldA,
        (mt::blas_scalar_t*)RS_val,
        (st::blas_scalar_t*)U_val,&ldU,
        (st::blas_scalar_t*)Vt_val,&ldVt,
        (st::blas_scalar_t*)&work,&lwork,
        (mt::blas_scalar_t*)rwork,iflag);
#else
    PHIST_TG_PREFIX(GESVD)((phist_blas_char*)(&jobu),(phist_blas_char*)(&jobvt),&m,&n,A_val,&ldA,RS_val,U_val,&ldU,Vt_val,&ldVt,work,&lwork,iflag);
#endif
    delete [] work;
    int ldS=svals_only?0: ldSigma;
    // copy the returned singular values to the diagonal of Sigma
    for (int i=0; i<mn; i++)
    {
      S_val[i*ldS+i]=(_ST_)(RS_val[i]);
    }
  }
}
