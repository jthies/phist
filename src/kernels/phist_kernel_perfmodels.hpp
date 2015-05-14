/*! \file phist_kernel_perfmodels.hpp
 * \brief defines performance models for kernel functions based on PHIST_PERFCHECK_VERIFY 
 *
*/
#ifndef PHIST_KERNEL_PERFMODELS_HPP
#define PHIST_KERNEL_PERFMODELS_HPP

#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_perfcheck.hpp"
#include "phist_kernels.h"



#ifdef PHIST_PERFCHECK

// define benchmarks
PHIST_PERFCHECK_BENCHMARK(STREAM_LOAD, phist_bench_stream_load);
PHIST_PERFCHECK_BENCHMARK(STREAM_TRIAD, phist_bench_stream_triad);
PHIST_PERFCHECK_BENCHMARK(STREAM_STORE, phist_bench_stream_store);
//PHIST_PERFCHECK_BENCHMARK(STREAM_FROM_DEVICE, phist_bench_stream_from_device);
//PHIST_PERFCHECK_BENCHMARK(STREAM_TO_DEVICE, phist_bench_stream_to_device);


//! checks performance for all functions that should require only neglectable time
#define PHIST_PERFCHECK_VERIFY_SMALL \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,0,0,0,0, 1.e-5);


//! checks performance of mvec_create
#define PHIST_PERFCHECK_VERIFY_MVEC_CREATE(map,nvec,iflag) \
  lidx_t n; \
  int nV = nvec; \
  PHIST_CHK_IERR(phist_map_get_local_length(map,&n,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_STORE(nV*n*sizeof(_ST_)));


//! checks performance of mvec_from_device
#define PHIST_PERFCHECK_VERIFY_MVEC_FROM_DEVICE(V,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_FROM_DEVICE(nV*n*sizeof(_ST_)));


//! checks performance of mvec_to_device
#define PHIST_PERFCHECK_VERIFY_MVEC_TO_DEVICE(V,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_TO_DEVICE(nV*n*sizeof(_ST_)));


//! checks performance of mvec_get_block
#define PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag) \
  lidx_t n; \
  int nV = jmax-jmin+1; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_TRIAD(2*nV*n*sizeof(_ST_)));


//! checks performance of mvec_set_block
#define PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag) \
  PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag)


//! checks performance of mvec_put_value
#define PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_STORE(nV*n*sizeof(_ST_)));


//! checks performance of mvec_dot_mvec (and mvec_norm2)
#define PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(V,W,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(V!=W),nV,0,0, STREAM_LOAD((nV+(V!=W)*nV)*n*sizeof(_ST_)));


//! checks performance of mvec_normalize
#define PHIST_PERFCHECK_VERIFY_MVEC_NORMALIZE(V,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_LOAD(nV*n*sizeof(_ST_))+STREAM_TRIAD(2*nV*n*sizeof(_ST_)));


//! checks performance of mvec_scale (and mvec_vscale)
#define PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0, STREAM_TRIAD(2*nV*n*sizeof(_ST_)));


//! checks performance of mvec_add_mvec
#define PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(a,X,b,Y,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(X,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nV,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=0),(b!=0),nV,0, STREAM_TRIAD(((a!=0)+(b!=0)+1)*nV*n*sizeof(_ST_)));


//! checks performance of mvec_vadd_mvec
#define PHIST_PERFCHECK_VERIFY_MVEC_VADD_MVEC(as,X,b,Y,iflag) \
  lidx_t n; \
  int nV; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(X,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nV,iflag),*iflag); \
  _MT_ a=mt::zero(); for(int i = 0; i < nV; i++) a=std::max(st::abs(as[i]),a); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=0),(b!=0),nV,0, STREAM_TRIAD(((a!=0)*nV+(b!=0)*nV+nV)*n*sizeof(_ST_)));


//! checks performance of mvecT_times_mvec
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag) \
  lidx_t n; \
  int nV, nW; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,0,0, STREAM_LOAD((nV+nW)*n*sizeof(_ST_)));


//! checks performance of mvec_times_sdMat
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(a,V,b,W,iflag) \
  lidx_t n; \
  int nV, nW; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=0),(b!=0),nV,nW, STREAM_TRIAD(((a!=0)*nV+(b!=0)*nW+nW)*n*sizeof(_ST_)));


//! checks performance of mvec_times_sdMat_inplace
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag) \
  lidx_t n; \
  int nV, nW; \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&nW,iflag),*iflag); \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,0,0, STREAM_TRIAD((nV+nW)*n*sizeof(_ST_)));


#else /* PHIST_PERFCHECK */


#define PHIST_PERFCHECK_VERIFY_SMALL
#define PHIST_PERFCHECK_VERIFY_MVEC_CREATE(map,nvec,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_FROM_DEVICE(V,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_TO_DEVICE(V,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(V,W,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_NORMALIZE(V,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(a,X,b,Y,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_VADD_MVEC(a,X,b,Y,iflag)
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(a,V,b,W,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag)


#endif


#endif /* PHIST_KERNEL_PERFMODELS_HPP */
