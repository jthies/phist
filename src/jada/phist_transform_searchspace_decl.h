//! apply a transformation matrix M to a given search space V and AV, BV and the projection H=V'AV
void SUBR(transform_searchSpace)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) M, bool generalizedEigenproblem, int* ierr);

//! in order to transform (e.g. shrink) a subspace we need a fast in-place operation for mvec_times_sdMat
void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, TYPE(sdMat_ptr) M_, lidx_t chunkSize, int* ierr);

