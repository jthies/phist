/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
// SVQB algorithm (Stathopoulos & Wu, SISC 23 (6),2165-2182).
// If the input vector has m columns and rank r, iflag=m-r is
// returned. Columns 0:iflag-1 of V will be zero, the remaining
// columns will be orthogonal. If the output matrix is denoted
// by Q, Q and V are related by Q=V*B. The third argument, E, 
// should be preallocated by the user with m elements, m being
// the number of columns in V. On successful return (iflag>=0),
// e[j] indicates the norm of V(:,j) before the orthogonali-  
// zation step. On exit, V(:,j), j>*iflag has 2-norm 1.
void SUBR(svqb)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) B, _MT_* D, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  bool robust = *iflag & PHIST_ROBUST_REDUCTIONS;
  *iflag=0;
  int m, rank;
  PHIST_CHK_IERR(SUBr(mvec_num_vectors)(V,&m,iflag),*iflag);
  _MT_ rankTol=mt::rankTol(robust);

  // S=V'V
  if( robust ) *iflag = PHIST_ROBUST_REDUCTIONS;
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),B,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_qb)(B,NULL,D,&rank,rankTol,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_to_device)(B,iflag),*iflag);
  if (robust) *iflag = PHIST_ROBUST_REDUCTIONS;
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,B,iflag),*iflag);
  *iflag=m-rank;
  return;
}
