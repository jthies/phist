/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

/*
#ifndef USE_FUSED_OPS
#define USE_FUSED_OPS 1
#endif
*/

//! orthogonalize an mvec against an already orthogonal one.
extern "C" void SUBR(orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_linearOp_ptr) B,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps,
                     int* rankVW,
                     int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  int iflag_in=*iflag;
#ifdef PHIST_HIGH_PRECISION_KERNELS_FORCE
  int robust=1;
#else
  int robust    =(*iflag&PHIST_ROBUST_REDUCTIONS);
#endif
  int randomize =(*iflag&PHIST_ORTHOG_RANDOMIZE_NULLSPACE);
  int m=0,k;
  // our current behavior is that if the user asks us to randomize the nullspace
  // explicitly, we try to determine the rank of the input [V, W]. Otherwise, we
  // just iterate until it has full rank. In JD, we do not explicitly ask for
  // randomization to avoid losing directions due to the 'rankTol'.
  _MT_ rankTol=randomize?mt::rankTol(robust):mt::zero();
  
  _MT_ orthoEps = mt::eps();
  if (robust) orthoEps*=std::sqrt(mt::eps());
  
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
  if (V!=NULL) PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);

    phist_const_comm_ptr comm=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(W,&comm,iflag),*iflag);
  
  TYPE(mvec_ptr) BW= W;
  
  TYPE(sdMat_ptr) WtW=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&WtW,k,k,comm,iflag),*iflag);
  phist::SdMatOwner<_ST_> _WtW(WtW);
  
  phist::MvecOwner< _ST_ > _BW;

  if (B!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_clone_shape)(&BW,W,iflag),*iflag);
    _BW.set(BW);
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    // note: we pass in WtW as the "VtW" argument because it should contain W'BW
    PHIST_CHK_IERR(B->fused_apply_mvTmv(st::one(),B->A,W,st::zero(),BW,NULL,WtW,iflag),*iflag);
  }
  else
  {
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),WtW,iflag),*iflag);
  }

  *iflag=iflag_in;
  SUBR(orthog_impl)(V,W,B,BW,WtW,R1,R2,numSweeps,rankVW,rankTol,orthoEps,iflag);
}

extern "C" void SUBR(orthog_impl)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_linearOp_ptr) B,
                     TYPE(mvec_ptr) BW,
                     TYPE(sdMat_ptr) WtW,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps,
                     int* rankVW,
                     _MT_ rankTol,
                     _MT_ orthoEps,
                     int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
#ifdef PHIST_HIGH_PRECISION_KERNELS_FORCE
  int robust=1;
#else
  int robust    =(*iflag&PHIST_ROBUST_REDUCTIONS);
#endif
  int randomize =(*iflag&PHIST_ORTHOG_RANDOMIZE_NULLSPACE);
  int triangular_R1=(*iflag&PHIST_ORTHOG_TRIANGULAR_R1);
  int dim0;

  int m=0,k;

  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
  if (V!=NULL) PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);

  if (V!=NULL)
  {
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
#ifdef USE_FUSED_OPS
    // the fused variant currently doesn't detect the rank of [V W] correctly in all cases, see issue #188
    PHIST_CHK_NEG_IERR(SUBR(orthogrrfused)(V, BW, R2, R1, WtW, iflag),*iflag);
#else
    PHIST_CHK_NEG_IERR(SUBR(orthogrrB)(V, W, BW, B, R2, R1, NULL,WtW, orthoEps,numSweeps,rankTol,iflag),*iflag);
#endif
    dim0=*iflag; // return value of orthog is rank of null space of [V W] on entry
    *rankVW=m+k-dim0;
  }
  else
  {
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
#ifdef USE_FUSED_OPS
    // fused orthog core doesn't allow V==NULL, so use orthogrr instead
    PHIST_CHK_NEG_IERR(SUBR(orthogrrfused)(V, W, BW, B, R2, R1, NULL,WtW, orthoEps,numSweeps,rankTol,iflag),*iflag);
#else
    PHIST_CHK_NEG_IERR(SUBR(orthogrrB)(V, W, BW, B, R2, R1, NULL,WtW, orthoEps,numSweeps,rankTol,iflag),*iflag);
#endif
    dim0=*iflag; // return value of orthog is rank of null space of [V W] on entry
    *rankVW=k-dim0;
  }
  int num_attempts=0;
  const int max_attempts=5;
  TYPE(sdMat_ptr) R1p=NULL,R1pp=NULL,R2p=NULL;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(SUBR(mvec_get_comm)(W,&comm,iflag),*iflag);
  if (randomize&&dim0>0)
  {
    if (m>0) PHIST_CHK_IERR(SUBR(sdMat_create)(&R2p,m,k,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1p,k,k,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1pp,k,k,comm,iflag),*iflag);
  }
  phist::SdMatOwner<_ST_> _R1pp(R1pp),_R1p(R1p), _R2p(R2p);
  while (randomize && dim0>0 && num_attempts++<max_attempts)
  {
    PHIST_SOUT(PHIST_DEBUG,"orthog: randomize %d column(s) to compensate rank deficiency (attempt %d)\n",
        dim0, num_attempts);
    // randomize last few columns, redo
    int rankW=k-dim0;
    TYPE(mvec_ptr) Wrnd=NULL;
    PHIST_CHK_IERR(SUBR(mvec_view_block)(W,&Wrnd,rankW,k-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_random)(Wrnd,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(Wrnd,iflag),*iflag);
  
    // compute W'W (or W'BW) with the new randomized column(s)
    if (B!=NULL)
    {
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(B->fused_apply_mvTmv(st::one(),B->A,W,st::zero(),BW,NULL,WtW,iflag),*iflag);
    }
    else
    {
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),WtW,iflag),*iflag);
    }

    // TODO: we do not use the fused core routine here either,
    //       it is not equipped for B-orthogonalization and I'm
    //       not sure about the consquences of issue #188 here,
    //       if the randomization fails and the rank is not    
    //       correctly determined, the result may be wrong.
    if (V!=NULL && false)
    {
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      SUBR(orthogrrfused)(V, BW, R2p, R1p, WtW, rankTol, iflag);
      dim0=*iflag;
    }
    else
    {
      // fused orthog core doesn't allow V==NULL, so use orthogrr instead
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      SUBR(orthogrrB)(V, W,BW,B,R2p, R1p, NULL,WtW, orthoEps,numSweeps,rankTol,iflag);
      dim0=*iflag; // return value of orthog is rank of null space of [V W] on entry
    }

    // update R1 and R2

    // we must not modify columns in R2 corresponding to random vectors!
    if (m>0)
    {
      TYPE(sdMat_ptr) R2p_r = NULL;
      PHIST_CHK_IERR(SUBR(sdMat_view_block)(R2p,&R2p_r,0,m-1,rankW,k-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R2p_r,st::zero(),iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p_r,iflag),*iflag);

      //R2=R2+R2'*R1;
      PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R2p,R1,st::one(),R2,iflag),*iflag);
    }
    
    //R1=R1p*R1;
    PHIST_CHK_IERR(SUBR(sdMat_add_sdMat)(st::one(),R1,st::zero(),R1pp,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_times_sdMat)(st::one(),R1p,R1pp,st::zero(),R1,iflag),*iflag);
    
    // zero-out the last dim0 columns of R2
    TYPE(sdMat_ptr) R1_r = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(R1,&R1_r,0,k-1,rankW,k-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R1_r,st::zero(),iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_r,iflag),*iflag);
  }
  
  // if the user explicitly asks for a triangular R1 factor, we have to permute
  // the final orthogonal Q (which is stored in W at this point).
  // We have Q*R1     = W - V*R2 = Q*P*R1p, so
  // first compute R1=P*R1p where R1p overwrites R1 and P is stored explicitly
  // then  update Q <- Q*P in-place.
  if (triangular_R1 && k>1)
  {
    if (R1p==NULL)
    {
      PHIST_CHK_IERR(SUBR(sdMat_create)(&R1p,k,k,comm,iflag),*iflag);
      _R1p.set(R1p);
    }
    TYPE(sdMat_ptr) P=R1p;
    PHIST_CHK_IERR(SUBR(sdMat_from_device)(R1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_qr)(P,R1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(R1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(P,iflag),*iflag);
    // update the output Q
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(W,P,iflag),*iflag);
    // for B-rothogonalization, also update BQ
    if (BW!=W) PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(BW,P,iflag),*iflag);
  }
  
  if (num_attempts==max_attempts)
  {
    *iflag=-8;
  }
  else if (*rankVW<m+k)
  {
    *iflag=+1;
  }
}
