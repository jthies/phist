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
  
  _MT_ orthoEps = mt::eps();
  
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
    PHIST_CHK_NEG_IERR(SUBR(orthogrrfused)(V, BW, R2, R1, WtW, iflag),*iflag);
    dim0=*iflag; // return value of orthog is rank of null space of W on entry
    *rankVW=m+k-dim0;
  }
  else
  {
    // fused orthog core doesn't allow V==NULL, so use orthogrr instead
    if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_NEG_IERR(SUBR(orthogrr)(V, BW, R2, R1, NULL,WtW, orthoEps,numSweeps,iflag),*iflag);
    dim0=*iflag; // return value of orthog is rank of null space of W on entry
    *rankVW=k-dim0;
  }

  int num_attempts=0;
  const int max_attempts=5;
  while (randomize && dim0>0 && num_attempts++<max_attempts)
  {
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

    if (V!=NULL)
    {
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      SUBR(orthogrrfused)(V, BW, R2, R1, WtW, iflag);
      dim0=m+k-*iflag;
    }
    else
    {
      // fused orthog core doesn't allow V==NULL, so use orthogrr instead
      if (robust) *iflag=PHIST_ROBUST_REDUCTIONS;
      SUBR(orthogrr)(V, BW, R2, R1, NULL,WtW, orthoEps,numSweeps,iflag);
      dim0=k-*iflag; // return value of orthog is rank of null space of W on entry
    }
    // zero-out the last dim0 columns of R2
    TYPE(sdMat_ptr) R1_r = NULL;
    PHIST_CHK_IERR(SUBR(sdMat_view_block)(R1,&R1_r,0,k-1,rankW,k-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R1_r,st::zero(),iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R1_r,iflag),*iflag);
  }
  
  if (B!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(st::one(),BW,st::zero(),W,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(BW,iflag),*iflag);
  }
  PHIST_CHK_IERR(SUBR(sdMat_delete)(WtW,iflag),*iflag);
  if (num_attempts==max_attempts)
  {
    *iflag=-8;
  }
  else if (*rankVW<m+k)
  {
    *iflag=+1;
  }
}



