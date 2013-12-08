//! in order to shrink a subspace we need a fast in-place operation for mvec_times_sdMat

void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, TYPE(sdMat_ptr) M_, lidx_t chunkSize, int* ierr);

