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
  int robust    =(*iflag&PHIST_ROBUST_REDUCTIONS);
  int randomize =(*iflag&PHIST_ORTHOG_RANDOMIZE_NULLSPACE);
  int m=0,k;
  
  _MT_ orthoEps = (_MT_)10.0*mt::eps();
  
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
  if (V!=NULL) PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);

    phist_const_map_ptr map=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(W,&map,iflag),*iflag);
    phist_const_comm_ptr comm=NULL;
    PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
  
  TYPE(mvec_ptr) BW= W;
  
  TYPE(sdMat_ptr) WtW=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_create)(&WtW,k,k,comm,iflag),*iflag);
  

  if (B!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_create)(&BW,map,k,iflag),*iflag);
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    // note: we pass in WtW as the "VtW" argument because it should contain W'BW
    PHIST_CHK_IERR(B->fused_apply_mvTmv(st::one(),B->A,W,st::zero(),BW,NULL,WtW,iflag),*iflag);
  }
  else
  {
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),WtW,iflag),*iflag);
  }
  
  int dim0;

  if (V!=NULL)
  {
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
// the fused variant currently doesn't detect the rank of [V W] correctly in all cases, see issue #188
//    PHIST_CHK_NEG_IERR(SUBR(orthogrrfused)(V, BW, R2, R1, WtW, iflag),*iflag);
    PHIST_CHK_NEG_IERR(SUBR(orthogrrB)(V, W, BW, B, R2, R1, NULL,WtW, orthoEps,numSweeps,iflag),*iflag);
    dim0=*iflag; // return value of orthog is rank of null space of [V W] on entry
    *rankVW=m+k-dim0;
  }
  else
  {
    // fused orthog core doesn't allow V==NULL, so use orthogrr instead
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_NEG_IERR(SUBR(orthogrrB)(V, W, BW, B, R2, R1, NULL,WtW, orthoEps,numSweeps,iflag),*iflag);
    dim0=*iflag; // return value of orthog is rank of null space of [V W] on entry
    *rankVW=k-dim0;
  }
  int num_attempts=0;
  const int max_attempts=5;
  TYPE(sdMat_ptr) R1p=NULL,R1pp=NULL,R2p=NULL;
  if (randomize&&dim0>0)
  {
    if (m>0) PHIST_CHK_IERR(SUBR(sdMat_create)(&R2p,m,k,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1p,k,k,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R1pp,k,k,comm,iflag),*iflag);
  }
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
      SUBR(orthogrrfused)(V, BW, R2p, R1p, WtW, iflag);
      dim0=*iflag;
    }
    else
    {
      // fused orthog core doesn't allow V==NULL, so use orthogrr instead
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      SUBR(orthogrrB)(V, W,BW,B,R2p, R1p, NULL,WtW, orthoEps,numSweeps,iflag);
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
  
  if (B!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(BW,iflag),*iflag);
  }
  PHIST_CHK_IERR(SUBR(sdMat_delete)(WtW,iflag),*iflag);
  if (R1p) PHIST_CHK_IERR(SUBR(sdMat_delete)(R1p,iflag),*iflag);
  if (R2p) PHIST_CHK_IERR(SUBR(sdMat_delete)(R2p,iflag),*iflag);
  if (num_attempts==max_attempts)
  {
    *iflag=-8;
  }
  else if (*rankVW<m+k)
  {
    *iflag=+1;
  }
}



