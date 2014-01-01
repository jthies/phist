//! apply a transformation matrix M to a given search space V and AV, BV and the projection H=V'AV
void SUBR(transform_searchSpace)(TYPE(mvec_ptr) V, TYPE(mvec_ptr) AV, TYPE(mvec_ptr) BV, TYPE(sdMat_ptr) H, TYPE(sdMat_ptr) M, bool generalizedEigenproblem, int *ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

  // get dimensions to create temporary storage
  lidx_t nV, minBase;
  const_comm_ptr_t comm;
  PHIST_CHK_IERR(SUBR( sdMat_get_nrows ) (M, &nV, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_get_ncols ) (M, &minBase, ierr), *ierr);

#ifdef PHIST_KERNEL_LIB_FORTRAN
  // temporary storage
  const_map_ptr_t map;
  PHIST_CHK_IERR( SUBR(mvec_get_map) (V, &map, ierr), *ierr);
  TYPE(mvec_ptr) Vtmp;
  PHIST_CHK_IERR( SUBR(mvec_create) (&Vtmp, map, minBase, ierr), *ierr);

  PHIST_CHK_IERR(SUBR( mvec_times_sdMat) (st::one(), V, M, st::zero(), Vtmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_set_block  ) (V, Vtmp, 0, minBase-1, ierr), *ierr);

  PHIST_CHK_IERR(SUBR( mvec_times_sdMat) (st::one(), AV, M, st::zero(), Vtmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_set_block  ) (AV, Vtmp, 0, minBase-1, ierr), *ierr);

  if( generalizedEigenproblem )
  {
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat) (st::one(), BV, M, st::zero(), Vtmp, ierr), *ierr);
    PHIST_CHK_IERR(SUBR( mvec_set_block  ) (BV, Vtmp, 0, minBase-1, ierr), *ierr);
  }

  PHIST_CHK_IERR( SUBR(mvec_delete)(Vtmp, ierr), *ierr);
#else
  // update V, AV and BV
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (V,  M, 64, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (AV, M, 64, ierr), *ierr);
  if( generalizedEigenproblem )
  {
    PHIST_CHK_IERR(SUBR( mvec_times_sdMat_inplace ) (BV, M, 64, ierr), *ierr);
  }
#endif

  // we need some communicator for ghost...
  PHIST_CHK_IERR(SUBR( mvec_get_comm ) (V, &comm, ierr), *ierr);
  TYPE(sdMat_ptr) H_ = NULL;
  TYPE(sdMat_ptr) Htmp = NULL;
  PHIST_CHK_IERR(SUBR( sdMat_create ) (&Htmp, nV, minBase, comm, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_view_block  )(H,    &H_,    0, minBase-1,    0, minBase-1,     ierr), *ierr);

  // update H <- M' * H * M
  PHIST_CHK_IERR(SUBR( sdMat_times_sdMat )(st::one(), H, M,    st::zero(), Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMatT_times_sdMat)(st::one(), M, Htmp, st::zero(), H_,   ierr), *ierr);

  // delete temp. storage
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (Htmp, ierr), *ierr);
  PHIST_CHK_IERR(SUBR( sdMat_delete ) (H_,   ierr), *ierr);
}


//! in order to shrink a subspace we need a fast in-place operation for mvec_times_sdMat
void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V_, TYPE(sdMat_ptr) M_, lidx_t chunkSize, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr = 0;

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
  ST *V = NULL;
  ST *M = NULL;
  lidx_t ldaV, ldaM;
  PHIST_CHK_IERR( SUBR( mvec_extract_view  ) (V_, &V,  &ldaV, ierr), *ierr);
  PHIST_CHK_IERR( SUBR( sdMat_extract_view ) (M_, &M,  &ldaM, ierr), *ierr);

#ifdef IS_COMPLEX
  typedef blas_cmplx_t BLAS_ST;
#else
  typedef ST BLAS_ST;
#endif

  int error = 0;
#pragma omp parallel reduction(||:error)
  {
#ifdef LIKWID_PERFMON
    likwid_markerStartRegion("mvec_times_sdMat_inplace");
#endif
    // create a buffer for a block of the mvec
    ST *work = (ST*)malloc(chunkSize*nvec*sizeof(ST));

    // blas needs int instead of lidx_t
    int i_chunkSize = chunkSize;
    int i_nvec = nvec;
    int i_ldaM = ldaM;
    int i_ldaV = ldaV;
    int info = 0;
    ST blas_ONE = st::one();
    ST blas_ZERO = st::zero();

    lidx_t i, j, off;
    lidx_t iChunk;
    lidx_t nChunks = nV/chunkSize;
#pragma omp for
    for(iChunk = 0; iChunk < nChunks; iChunk++)
    {
      off = iChunk*chunkSize;
      // copy block to work
      for(j = 0; j < nvec; j++)
        for(i = 0; i < chunkSize; i++)
          work[chunkSize*j+i] = V[ldaV*j+off+i];

      // call gemm with work buffer
      PREFIX(GEMM) ("N", "N", &i_chunkSize, &i_nvec, &i_nvec, (BLAS_ST*)&blas_ONE, (BLAS_ST*)work, &i_chunkSize, (BLAS_ST*)M, &i_ldaM, (BLAS_ST*)&blas_ZERO, (BLAS_ST*)&(V[off]), &i_ldaV, &info);
      error = error || (info != 0);
    }


#pragma omp single
    {
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
        PREFIX(GEMM) ("N", "N", &i_chunkSize, &i_nvec, &i_nvec, (BLAS_ST*)&blas_ONE, (BLAS_ST*)work, &i_chunkSize, (BLAS_ST*)M, &i_ldaM, (BLAS_ST*)&blas_ZERO, (BLAS_ST*)&(V[off]), &i_ldaV, &info);
        error = error || (info != 0);
      }
    }


    free(work);
#ifdef LIKWID_PERFMON
    likwid_markerStopRegion("mvec_times_sdMat_inplace");
#endif
  }

  PHIST_CHK_IERR(*ierr = error ? -1: 0, *ierr);
}

