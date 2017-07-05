/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

// implementation of gmres on several systems simultaneously
void SUBR(blockedQMR_iterate)(TYPE(const_linearOp_ptr) Aop, TYPE(const_linearOp_ptr) Pop,
        TYPE(const_mvec_ptr) rhs, TYPE(mvec_ptr) sol, int numSys, int* nIter, _MT_ tol[], int* iflag)
{
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  if (numSys==0) return; // do not appear in timing stats
  PHIST_ENTER_FCN(__FUNCTION__);

  int maxIter = (*nIter)>0 ? *nIter: 9999999;
  *nIter=0;

  *iflag=PHIST_NOT_IMPLEMENTED;

}
