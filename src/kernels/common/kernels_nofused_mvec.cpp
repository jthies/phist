void SUBR(mvec_times_sdMat_aug)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                  TYPE(const_sdMat_ptr) C,
                                      _ST_ beta,  TYPE(mvec_ptr)        W,
                                                  TYPE(sdMat_ptr)       D,
                                                  int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(mvec_times_sdMat)(alpha,V,C,beta,W,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

void SUBR(mvecT_times_mvec_times_sdMat_inplace)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                            TYPE(mvec_ptr)        W,
                                                            TYPE(const_sdMat_ptr) C,
                                                _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // call two kernels
  // don't freeze inside here if the first one returns an error on only some processes
  int iflag1 = *iflag;
  SUBR(mvec_times_sdMat_inplace)(W,C,&iflag1);
  int iflag2 = *iflag;
  SUBR(mvecT_times_mvec)(alpha,V,W,beta,D,&iflag2);
  PHIST_CHK_IERR(*iflag = iflag1,*iflag);
  PHIST_CHK_IERR(*iflag = iflag2,*iflag);
}

void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V, 
                                                 TYPE(const_sdMat_ptr) C,
                                                 TYPE(mvec_ptr) W, 
                                                 TYPE(const_sdMat_ptr) D,
                                                 int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  // remember flags
  int flags = *iflag;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(W,D,iflag),*iflag);
  *iflag = flags;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),V,C,st::one(),W,iflag),*iflag);
}


// ghost implements this one but not the above right now
#ifndef PHIST_KERNEL_LIB_GHOST

// augmented spMVM kernel available in GHOST

// like sparseMat_times_mvec_add_mvec, followed by z=a*y+b*z. if z!=NULL.
// if dot_xx!=NULL, it will contain mvec_dot_mvec(x,x) on output
// and similarly for dot_xy and dot_yy (final y being used)
//
// Kernel libraries that do not offer this can include 
// common/kernels_no_fused.cpp for a fallback variant
void SUBR(sparseMat_times_mvec_vaug)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ const shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, 
        _ST_ a, _ST_ b, TYPE(mvec_ptr) z,
        _ST_* dot_xx, _ST_* dotxy, _ST_* dotyy, 
        int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A,shifts,x,beta,y,iflag);
  if (z!=NULL)
  {
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(a,y,b,z,iflag),*iflag);
  }
  if (dot_xx)
  {
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(x,x,dot_xx,iflag);
  }
  if (dot_xy)
  {
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(x,y,dot_xy,iflag);
  }
  if (dot_yy)
  {
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(x,x,dot_yy,iflag);
  }
}
#endif
