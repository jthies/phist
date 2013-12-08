//! in order to shrink a subspace we need a fast in-place operation for mvec_times_sdMat

void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, TYPE(sdMat_ptr) M_, lidx_t chunkSize, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
#ifdef LIKWID_PERFMON
  likwid_markerStartRegion("mvec_times_sdMat_inplace");
#endif

  // get dimensions
  lidx_t nV, nvec, nM, mM;
  PHIST_CHK_IERR( SUBR( mvec_my_length     ) (V_, &nV,        ierr), *ierr);
  PHIST_CHK_IERR( SUBR( mvec_num_vectors   ) (V_, &nvec,      ierr), *ierr);
  PHIST_CHK_IERR( SUBR( sdMat_get_nrows    ) (M_, &nM,        ierr), *ierr);
  PHIST_CHK_IERR( SUBR( sdMat_get_ncols    ) (M_, &mM,        ierr), *ierr);

  // check dimensions
  PHIST_CHK_IERR(*ierr = (nvec == nM ? 0 : -1), *ierr);
  PHIST_CHK_IERR(*ierr = (nvec >= mM ? 0 : -1), *ierr);

  // get raw view on V and M
  _ST_ *V = NULL;
  _ST_ *M = NULL;
  lidx_t ldaV, ldaM;
  PHIST_CHK_IERR( SUBR( mvec_extract_view  ) (V_, &V,  &ldaV, ierr), *ierr);
  PHIST_CHK_IERR( SUBR( sdMat_extract_view ) (M_, &M,  &ldaM, ierr), *ierr);

  // create a buffer for a block of the mvec
  _ST_ *work = (_ST_*)malloc(chunkSize*nvec*sizeof(_ST_));

  // blas needs int instead of lidx_t
  int i_chunkSize = chunkSize;
  int i_nvec = nvec;
  int i_ldaM = ldaM;
  int i_ldaV = ldaV;
#ifdef IS_COMPLEX
  _ST_ blas_ONE_ = ONE;
  _ST_ blas_ZERO_ = ZERO;
  const _ST_* blas_ONE = &blas_ONE_;
  const _ST_* blas_ZERO = &blas_ZERO_;
#else
  _ST_ blas_ONE = ONE;
  _ST_ blas_ZERO = ZERO;
#endif

  lidx_t i, j, k, off;
  lidx_t iChunk;
  lidx_t nChunks = nV/chunkSize;
  for(iChunk = 0; iChunk < nChunks; iChunk++)
  {
    off = iChunk*chunkSize;
    // copy block to work
    for(j = 0; j < nvec; j++)
      for(i = 0; i < chunkSize; i++)
        work[chunkSize*j+i] = V[ldaV*j+off+i];

    // call gemm with work buffer
    PHIST_CHK_IERR( PREFIX(GEMM) ("N", "N", &i_chunkSize, &i_nvec, &i_nvec, &blas_ONE, work, &i_chunkSize, M, &i_ldaM, &blas_ZERO, &(V[off]), &i_ldaV, ierr), *ierr);
  }


  // do the remaining part
  off = nChunks*chunkSize;
  chunkSize = nV-nChunks*chunkSize;
  i_chunkSize = chunkSize;

  if( chunkSize > 0 )
  {
    // copy block to work
    for(j = 0; j < nvec; j++)
      for(i = 0; i < chunkSize; i++)
          work[chunkSize*j+i] = V[ldaV*j+off+i];

    // call gemm with work buffer
    PHIST_CHK_IERR( PREFIX(GEMM) ("N", "N", &i_chunkSize, &i_nvec, &i_nvec, &blas_ONE, work, &i_chunkSize, M, &i_ldaM, &blas_ZERO, &(V[off]), &i_ldaV, ierr), *ierr);
  }


  free(work);
#ifdef LIKWID_PERFMON
  likwid_markerStopRegion("mvec_times_sdMat_inplace");
#endif
}

