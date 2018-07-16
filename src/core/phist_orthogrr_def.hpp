/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

// we can use either a Cholesky factorization of the Gram matrix V'BV (CholQR) or the SVQB algorithm,
// the required transformation matrices are defined in the following two helpr functions.
//
// Our orthog routines (e.g. orthogrr below) should deliver the relation
//
//      Q*R_1 = V - W*R_2 with R_2=W'*V
//
// However, we can only construct R_1 as the inverse of the Cholesky factor, so we use CholQR here.
// The core of the SVQB algorithm delivers a matrix B s.t. Q=V*B, and the left pseudo-inverse B+ s.t.
// (B+)B=[I 0; 0 0] with zeros if V does not have full rank (a \rank-r identity matrix'), so we can't
// recover V from Q=VB and B+.

// On input: RR=V'V (the m x m Gramian 'M')
// returns: RR and R_1 s.t. M=RR'*RR, R_1*RR is a 'rank-r identity matrix', and Q=V*R_1 is orthonormal
// (up to the last m-r columns, r being the rank of V on input), which are zero.
static void SUBR(orthogrr_cholrr)(TYPE(sdMat_ptr) RR, TYPE(sdMat_ptr) R_1, int* rank, _MT_ rankTol, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  int m = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(RR,&m,iflag),*iflag);
  // create R and R_1 from RR
  *rank = 0;
  {
    // stable rank-revealing cholesky
    int perm[m];
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(RR,perm,rank,rankTol,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(RR,perm,*rank,R_1,iflag),*iflag);
  }
//PHIST_CHK_IERR(SUBR(sdMat_print)(R_1,iflag),*iflag);
}

// given the Gramian M=V'Vv in RR, returns B s.t. Q=V*B is orthonormal (up to rank deficiency).
// M will be overwritten with the right pseudo-inverse of B s.t. M*B on exit is a 'rank r identity matrix'.
// Note that we cannot reproduce V as a product of Q and the pinv of B in the rank-deficient case.
static void SUBR(orthogrr_svqb)(TYPE(sdMat_ptr) RR, TYPE(sdMat_ptr) R_1, int* rank, _MT_ rankTol, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // compute B s.t. V*B=Q is orthogonal, with RR=V'V on input
  // and pinv(B) on output, and R_1=B on output.
  
  // we first copy the input matrix because B and B_1 are exchanged in the definition of
  // the kernel routine sdMat_qb:
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),RR,st::zero(),R_1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_qb)(R_1,RR,NULL,rank,rankTol,iflag),*iflag);
}


static void SUBR(sdMat_rank_identity)(TYPE(sdMat_ptr) I, const int k, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(I,st::zero(),iflag),*iflag);
  TYPE(sdMat_ptr) Ik = NULL;
  if (k>0)
  {
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(I,&Ik,0,k-1,0,k-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_identity)(Ik,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(Ik,iflag),*iflag);
  }
}

//
void SUBR(orthogrr)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, TYPE(const_sdMat_ptr) WtW_I, TYPE(sdMat_ptr) VtV, 
_MT_ desiredEps, int maxIter, _MT_ rankTol, int* iflag)
{
  SUBR(orthogrrB)(W,V,V,NULL,R2,R1,WtW_I,VtV,desiredEps,maxIter,rankTol,iflag);
}

// Q*R_1 = V - W*R_2 with R_2=W'*V
// result returned in place as V
// correct V'V must be supplied and is returned!
void SUBR(orthogrrB)(TYPE(const_mvec_ptr) W, TYPE(mvec_ptr) V, TYPE(mvec_ptr) BV, TYPE(const_linearOp_ptr) B_op,
        TYPE(sdMat_ptr) R2, TYPE(sdMat_ptr) R1, TYPE(const_sdMat_ptr) WtW_I, TYPE(sdMat_ptr) VtV, \
        _MT_ desiredEps, int maxIter, _MT_ rankTol, int* iflag)
{
#include "phist_std_typedefs.hpp"
    PHIST_ENTER_FCN(__FUNCTION__);
    bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
    *iflag=0;
    // get dimensions
    int m, k = 0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    if( W != NULL ) {PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);}
    phist_const_comm_ptr comm;
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
    // sanity checks - if your B operator is the identity matrix, just leave it away and set W==BW
    PHIST_CHK_IERR(*iflag=(B_op!=NULL && V==BV)?PHIST_INVALID_INPUT:0,*iflag);
    PHIST_CHK_IERR(*iflag=(B_op==NULL && V!=BV)?PHIST_INVALID_INPUT:0,*iflag);

    if( k > 0 && WtW_I )
    {
      PHIST_CHK_IERR(SUBR(sdMat_create)(&WtW_inv,k,k,comm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_create)(&EWtV,k,m,comm,iflag),*iflag);
      // construct cholesky factorization of WtW
      PHIST_CHK_IERR(SUBR(sdMat_identity)(WtW_inv,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),WtW_I,st::one(),WtW_inv,iflag),*iflag);
      permWtW = new int[k];
      rankWtW = 0;
      PHIST_CHK_IERR(SUBR(sdMat_cholesky)(WtW_inv,permWtW,&rankWtW,rankTol,iflag),*iflag);
    }

    // memory management
    SdMatOwner<_ST_> _R(R),_R1_tmp(R1_tmp),_R_1(R_1),_WtV(WtV),_VtV_I(VtV_I),_WtW_inv(WtW_inv),_EWtV(EWtV);

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
    bool VtV_updated = true;
    for(; iter < maxIter; iter++)
    {
      // we already have the current VtV!
      // so calculate first "R" factor
      int Vrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
*iflag=robust?PHIST_ROBUST_REDUCTIONS:0;
#ifdef ORTHOGRR_USE_SVQB
      PHIST_CHK_IERR(SUBR(orthogrr_svqb)(R,R_1,&Vrank,rankTol,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&Vrank,rankTol,iflag),*iflag);
#endif
      rank = std::min(rank,Vrank);
      if( k == 0 )
      {
        iter++;
        break;
      }

      // BV <- BV*R_1 and WtV <- W'*BV
      PHIST_CHK_IERR(SUBR(sdMat_to_device)(R_1,iflag),*iflag);
      if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(fused_mvsdi_mvTmv)(st::one(),W,BV,R_1,st::zero(),WtV,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(WtV,iflag),*iflag);
      if (B_op!=NULL)
      {
        if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,R_1,iflag),*iflag);
      }
      if( WtW_I )
      {
        // correct WtV by error of WtW from
        // 0 = W^T(V-W(WtV+X))
        //   = WtV - WtW (WtV+X)
        //   = WtV - (I+E) (WtV+X)
        //   = -E WtV - (I+E)X
        // So we obtain the correction
        // X = - WtW^(-1) WtW_I WtV
        PHIST_CHK_IERR(SUBR(sdMat_from_device)(WtW_inv,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtW_I,WtV,st::zero(),EWtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(WtW_inv,permWtW,rankWtW,EWtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_forwardSubst_sdMat)(WtW_inv,permWtW,rankWtW,EWtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(-st::one(),EWtV,st::one(),WtV,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_to_device)(WtV,iflag),*iflag);
      }

      // update R1 <- R * R1
      if( R1 != NULL )
      {
        PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1_tmp,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R1_tmp,st::zero(),R1,iflag),*iflag);
      }
      // check if we really need another orthogonalization step with W
      if( robust ) 
      {
        VtV_err = std::max(VtV_err*mt::eps(),mt::eps());
      }
      else
      {
        VtV_err = std::max(VtV_err*mt::sqrt(mt::eps()),mt::eps());
      }
      PHIST_CHK_IERR(SUBR(sdMat_normF)(WtV,&WtV_err,iflag),*iflag);
      VtV_err += WtV_err;
      PHIST_SOUT(PHIST_EXTREME, "orthogRR: iter %d phase 1, desired eps %8.4e, WtV err. %8.4e, (est.) VtV err. %8.4e\n", iter, desiredEps, WtV_err, VtV_err);
      if( WtV_err <= desiredEps && WtV_err*VtV_err <= desiredEps*desiredEps )
      {
        // calculated WtV_err is small enough (in contrast to estimated WtV_err from the previous iteration)
        // we may need to recalculate VtV below
        VtV_updated = false;
        break;
      }

      // update R2 <- R2 + W'*V = R2 + WtV*R
      if( R2 != NULL ) {PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),WtV,R,st::one(),R2,iflag),*iflag);}

      // V <- V - W*WtV, updating VtV. TODO: fused kernel with B
      if (B_op==NULL)
      {
        if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(SUBR(fused_mvsd_mvTmv)(-st::one(),W,WtV,st::one(),V,VtV,iflag),*iflag);
      }
      else
      {
        if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),W,WtV,st::one(),V,iflag),*iflag);
        if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
        PHIST_CHK_IERR(B_op->fused_apply_mvTmv(st::one(),B_op->A,V,st::zero(),BV,NULL,VtV,iflag),*iflag);
      }

      // calculate new R factor
      int WVrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(VtV,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
      *iflag=robust?PHIST_ROBUST_REDUCTIONS:0;
#ifdef ORTHOGRR_USE_SVQB
      PHIST_CHK_IERR(SUBR(orthogrr_svqb)(R,R_1,&WVrank,rankTol,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&WVrank,rankTol,iflag),*iflag);
#endif
      rank = std::min(rank,WVrank);

      // estimate error
      PHIST_CHK_IERR(SUBR(sdMat_normF)(R_1,&WtV_err,iflag),*iflag);
      WtV_err *= mt::eps();
      PHIST_CHK_IERR(SUBR(sdMat_rank_identity)(VtV_I,rank,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,-st::one(),VtV_I,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_normF)(VtV_I,&VtV_err,iflag),*iflag);

      // possibly perform another full sweep with W, if R_1 is too badly conditioned.
      PHIST_SOUT(PHIST_EXTREME, "orthogRR: iter %d phase 2, desired eps %8.4e, (est.) WtV err. %8.4e, VtV err. %8.4e\n", iter, desiredEps, WtV_err, VtV_err);
      if( WtV_err <= desiredEps )
      {
        iter++;
        break;
      }
    }

    // continue iteration with VtV
    for(; iter <= maxIter; iter++)
    {
      if( VtV_err <= desiredEps ) break;

      // we need to recalculate VtV
      if( !VtV_updated )
      {
        if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
        if (B_op!=NULL)
        {
          PHIST_CHK_IERR(B_op->fused_apply_mvTmv(st::one(),B_op->A,V,st::zero(),BV,NULL,VtV,iflag),*iflag);
        }
        else
        {
          PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),VtV,iflag),*iflag);
        }
        VtV_updated = true;
        PHIST_CHK_IERR(SUBR(sdMat_from_device)(VtV,iflag),*iflag);
      }
      else
      {
        // V <- V*R_1, VtV <- V'*BV. TODO: fused kernel with B
        PHIST_CHK_IERR(SUBR(sdMat_to_device)(R_1,iflag),*iflag);
        if (B_op!=NULL)
        {
          if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
          PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(BV,R_1,iflag),*iflag);
          if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
          PHIST_CHK_IERR(SUBR(fused_mvsdi_mvTmv)(st::one(),BV,V,R_1,st::zero(),VtV,iflag),*iflag);
          
        }
        else
        {
          if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
          PHIST_CHK_IERR(SUBR(fused_mvsdi_mvTmv)(st::one(),V,V,R_1,st::zero(),VtV,iflag),*iflag);
        }

        // update R1 <- R * R1
        if( R1 != NULL )
        {
          PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1_tmp,iflag),*iflag);
          PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R,R1_tmp,st::zero(),R1,iflag),*iflag);
        }
        PHIST_CHK_IERR(SUBR(sdMat_from_device)(VtV,iflag),*iflag);
      }

      // calculate new R factor
      int Vrank = 0;
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,st::zero(),R,iflag),*iflag);
      *iflag=robust?PHIST_ROBUST_REDUCTIONS:0;
#ifdef ORTHOGRR_USE_SVQB
      PHIST_CHK_IERR(SUBR(orthogrr_svqb)(R,R_1,&Vrank,rankTol,iflag),*iflag);
#else
      PHIST_CHK_IERR(SUBR(orthogrr_cholrr)(R,R_1,&Vrank,rankTol,iflag),*iflag);
#endif
      rank = std::min(rank,Vrank);

      // calculate error
      PHIST_CHK_IERR(SUBR(sdMat_rank_identity)(VtV_I,rank,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),VtV,-st::one(),VtV_I,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_normF)(VtV_I,&VtV_err,iflag),*iflag);
      PHIST_SOUT(PHIST_EXTREME, "orthogRR: iter %d phase 3, desired eps %8.4e, VtV err. %8.4e\n", iter-1, desiredEps, VtV_err);
    }

    // upload resulting sdMats
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(R1,iflag),*iflag);
    if (R2!=NULL) PHIST_CHK_IERR(SUBR(sdMat_to_device)(R2,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(VtV,iflag),*iflag);
    

    // the return value of this function is the rank of the null space of V on entry
    *iflag=m-rank;
}
