#include "phist_config.h"
#include "phist_void_aliases.h"
#include "phist_kernels.h"

extern "C" void SUBR(sparseMat_get_context)(TYPE(const_sparseMat_ptr) A, phist_const_context_ptr *vctx, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  phist_default_context* ctx=contextCollection[A];
  if (ctx!=NULL) return; // context has been created before -> return it again
  ctx=new phist_default_context;
  PHIST_CHK_IERR(SUBR(sparseMat_get_row_map)(A,&ctx->row_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&ctx->domain_map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&ctx->range_map,iflag),*iflag);
  contextCollection[A]=ctx;
  *vctx=(phist_const_context_ptr)ctx;
}
