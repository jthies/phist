/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
#include "phist_void_aliases.h"
#include "phist_kernels.h"
#include "./default_context.h"

extern "C" void SUBR(sparseMat_get_context)(TYPE(const_sparseMat_ptr) A, phist_const_context_ptr *vctx, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  phist::internal::default_context* ctx=phist::internal::get_default_context(A);
  PHIST_CHK_IERR(SUBR(sparseMat_get_row_map)(A,&ctx->row_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_col_map)(A,&ctx->col_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&ctx->domain_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&ctx->range_map,iflag),*iflag);
  *vctx=(phist_const_context_ptr)ctx;
}
