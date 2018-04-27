/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file phist_kernel_perfmodels.hpp
 * \brief defines performance models for kernel functions based on PHIST_PERFCHECK_VERIFY 
 *
*/
#ifndef PHIST_KERNEL_PERFMODELS_HPP
#define PHIST_KERNEL_PERFMODELS_HPP

#include "phist_config.h"
#ifndef DOXYGEN
#include "phist_tools.h"
#include "phist_perfcheck.hpp"
#include "phist_kernels.h"
#endif


#ifdef PHIST_PERFCHECK

#ifdef PHIST_PERFCHECK_MAP_LENGTH
#undef PHIST_PERFCHECK_MAP_LENGTH
#endif
#ifdef PHIST_PERFCHECK_MVEC_LENGTH
#undef PHIST_PERFCHECK_MVEC_LENGTH
#endif
#ifdef PHIST_PERFCHECK_MVEC_LEN_T
#undef PHIST_PERFCHECK_MVEC_LEN_T
#endif

#ifdef PHIST_PERFCHECK_SEPARATE_OUTPUT
#define PHIST_PERFCHECK_MVEC_LENGTH SUBR(mvec_my_length)
#define PHIST_PERFCHECK_MAP_LENGTH phist_map_get_local_length
#define PHIST_PERFCHECK_MVEC_LEN_T phist_lidx
#else
#define PHIST_PERFCHECK_MVEC_LENGTH SUBR(mvec_global_length)
#define PHIST_PERFCHECK_MAP_LENGTH phist_map_get_global_length
#define PHIST_PERFCHECK_MVEC_LEN_T phist_gidx
#endif

// flops for operation Y <- a*X + b*Y, including some special cases like a=1, b=0 etc.
#ifndef PHIST_PERFCHECK_AXPBY_FLOPS
#define PHIST_PERFCHECK_AXPBY_FLOPS(N,a,b)  (N*((a!=_ST_(0))*(1+(a!=_ST_(1)))+(b!=_ST_(0))*(1+(b!=_ST_(1)))))
#endif

// define benchmarks
PHIST_PERFCHECK_BENCHMARK(STREAM_LOAD, phist_bench_stream_load);
PHIST_PERFCHECK_BENCHMARK(STREAM_TRIAD, phist_bench_stream_triad);
PHIST_PERFCHECK_BENCHMARK(STREAM_STORE, phist_bench_stream_store);
//PHIST_PERFCHECK_BENCHMARK(STREAM_FROM_DEVICE, phist_bench_stream_from_device);
//PHIST_PERFCHECK_BENCHMARK(STREAM_TO_DEVICE, phist_bench_stream_to_device);
#define NO_PERFMODEL_AVAILABLE_YET 0.0


//! checks performance for all functions that should require only neglectable time
#define PHIST_PERFCHECK_VERIFY_SMALL \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,0,0,0,0,0,0,0, 1.e-5,0);


//! prints a warning that the performance model is not yet implemented
#define PHIST_PERFCHECK_NOT_IMPLEMENTED(flops) \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,0,0,0,0,0,0,0,NO_PERFMODEL_AVAILABLE_YET,flops);


//! performance model for sparse matrix-vector product with fused operations:
//! Y <- beta*Y + (A + shift*I)X, s<-gamma*X'Y, t<-delta*Y'Y
#define PHIST_PERFCHECK_VERIFY_SPMV(a,A,shift,X,b,Y,g,d,num_dots,iflag) \
  int _nV; \
  int64_t _nnz; \
  phist_lidx _n; \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&_nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_my_length)(X,&_n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(sparseMat_get_local_nnz)(A,&_nnz,iflag),*iflag); \
  double flops = (a!=0)*double(2*_nnz*_nV + (shift!=_ST_(0)*_n*_nV)) + (b!=_ST_(0))*double((1+(b!=_ST_(1)))*_n*_nV) + PHIST_PERFCHECK_AXPBY_FLOPS(_n*_nV,g,d) + double(2*num_dots*_nV*_n); \
PHIST_PERFCHECK_NOT_IMPLEMENTED(flops);

//! checks performance of mvec_create
#define PHIST_PERFCHECK_VERIFY_MVEC_CREATE(map,nvec,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV = nvec; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MAP_LENGTH(map,&n,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_STORE(nV*n*sizeof(_ST_)),0);


//! checks performance of mvec_from_device
#define PHIST_PERFCHECK_VERIFY_MVEC_FROM_DEVICE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_FROM_DEVICE(nV*n*sizeof(_ST_)),0);


//! checks performance of mvec_to_device
#define PHIST_PERFCHECK_VERIFY_MVEC_TO_DEVICE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_TO_DEVICE(nV*n*sizeof(_ST_)),0);


//! \def PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag)
//! checks performance of mvec_get_block
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV = jmax-jmin+1; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_TRIAD(2*nV*n*sizeof(_ST_)),0);

#else

// realistic model which respects cache line length
// note: do not check iflag on mvec_extract_view because it may
// fail but set lda correctly for ghost+cuda
#define PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV = jmax-jmin+1; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  _ST_ *V_raw, *Vb_raw; \
  phist_lidx ldV, ldVb; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))Vblock,&Vb_raw,&ldVb,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nVb_ = std::min(ldVb,((nV-1)/cl_size+1)*cl_size); \
  if( nVb_+cl_size > ldVb ) nVb_ = ldVb; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,ldV,ldVb,(nV!=nVb_),0,0,0,STREAM_TRIAD((nV_+nVb_+(nV!=nVb_)*nVb_)*n*sizeof(_ST_)),0);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag)
//! checks performance of mvec_set_block
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag) \
  PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag)

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV = jmax-jmin+1; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  _ST_ *V_raw, *Vb_raw; \
  phist_lidx ldV, ldVb; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))Vblock,&Vb_raw,&ldVb,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nVb_ = std::min(ldVb,((nV-1)/cl_size+1)*cl_size); \
  if( nVb_+cl_size > ldVb ) nVb_ = ldVb; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,ldV,ldVb,(nV!=nV_),0,0,0,STREAM_TRIAD((nV_+nVb_+(nV!=nV_)*nV_)*n*sizeof(_ST_)),0);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag)
//! checks performance of mvec_put_value
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_STORE(nV*n*sizeof(_ST_)),0);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  _ST_ *V_raw; \
  phist_lidx ldV; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,ldV,(nV!=nV_),0,0,0,0, (nV==nV_)*STREAM_STORE(nV*n*sizeof(_ST_)) + (nV!=nV_)*STREAM_TRIAD(2*nV_*n*sizeof(_ST_)),0);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(V,W,iflag)
//! checks performance of mvec_dot_mvec (and mvec_norm2)
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(V,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(V!=W),nV,0,0,0,0,0, STREAM_LOAD((nV+(V!=W)*nV)*n*sizeof(_ST_)),2*n);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(V,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  _ST_ *V_raw, *W_raw; \
  phist_lidx ldV, ldW; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))W,&W_raw,&ldW,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nW_ = std::min(ldW,((nV-1)/cl_size+1)*cl_size); \
  if( nW_+cl_size > ldW ) nW_ = ldW; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(V!=W),nV,ldV,ldW,0,0,0, STREAM_LOAD((nV_+(V!=W)*nW_)*n*sizeof(_ST_)),2*n);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_NORMALIZE(V,iflag)
//! checks performance of mvec_normalize
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_NORMALIZE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_LOAD(nV*n*sizeof(_ST_))+STREAM_TRIAD(2*nV*n*sizeof(_ST_)),3*n);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_NORMALIZE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  _ST_ *V_raw; \
  phist_lidx ldV; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,ldV,0,0,0,0,0, STREAM_LOAD(nV_*n*sizeof(_ST_))+STREAM_TRIAD(2*nV_*n*sizeof(_ST_)),3*n);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag)
//! checks performance of mvec_scale (and mvec_vscale)
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,0,0,0,0,0, STREAM_TRIAD(2*nV*n*sizeof(_ST_)),n);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  _ST_ *V_raw; \
  phist_lidx ldV; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,ldV,0,0,0,0,0, STREAM_TRIAD(2*nV_*n*sizeof(_ST_)),n);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(a,X,b,Y,iflag)
//! checks performance of mvec_add_mvec
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(a,X,b,Y,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(X,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nV,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=_ST_(0)),(b!=_ST_(0)),nV,0,0,0,0, STREAM_TRIAD(((a!=_ST_(0))+(b!=_ST_(0))+1)*nV*n*sizeof(_ST_)), PHIST_PERFCHECK_AXPBY_FLOPS(n*nV,a,b));

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(a,X,b,Y,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(X,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nV,iflag),*iflag); \
  _ST_ *X_raw, *Y_raw; \
  phist_lidx ldX, ldY; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))X,&X_raw,&ldX,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))Y,&Y_raw,&ldY,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nX_ = std::min(ldX,((nV-1)/cl_size+1)*cl_size); \
  if( nX_+cl_size > ldX ) nX_ = ldX; \
  int nY_ = std::min(ldY,((nV-1)/cl_size+1)*cl_size); \
  if( nY_+cl_size > ldY ) nY_ = ldY; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=_ST_(0)),(b!=_ST_(0)),nV,nX_,nY_,0,0, STREAM_TRIAD(((a!=_ST_(0))*nX_+(b!=_ST_(0))*nY_+nY_)*n*sizeof(_ST_)),PHIST_PERFCHECK_AXPBY_FLOPS(n*nV,a,b));

#endif


//! checks performance of mvec_vadd_mvec
#define PHIST_PERFCHECK_VERIFY_MVEC_VADD_MVEC(as,X,b,Y,iflag) \
  int tmp_iflag_ = *iflag; \
  _MT_ as_=mt::zero(); \
  { \
    int nV; \
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X,&nV,iflag),*iflag); \
    for(int i = 0; i < nV; i++) as_=std::max(st::abs(as[i]),as_); \
  } \
  _ST_ a_ = as_; \
  *iflag = tmp_iflag_; \
  PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(a_,X,b,Y,iflag);


//! \def PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag)
//! checks performance of mvecT_times_mvec
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,(V!=W),*iflag,0,0,0, STREAM_LOAD((nV+(V!=W)*nW)*n*sizeof(_ST_)),2*n*nV*nV);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  _ST_ *V_raw, *W_raw; \
  phist_lidx ldV, ldW; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))W,&W_raw,&ldW,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nW_ = std::min(ldW,((nW-1)/cl_size+1)*cl_size); \
  if( nW_+cl_size > ldW ) nW_ = ldW; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,(V!=W),ldV,ldW,*iflag,0, STREAM_LOAD((nV_+(V!=W)*nW_)*n*sizeof(_ST_)),2*n*nV*nV);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC_TIMES_SDMAT(V,W,iflag)
//! checks performance of fused_mvsdi_mvTmv
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC_TIMES_SDMAT(V,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,(V!=W),*iflag,0,0,0, STREAM_TRIAD((nV+(V!=W)*nW+nW)*n*sizeof(_ST_)),4*n*nV*nV);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC_TIMES_SDMAT(V,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  _ST_ *V_raw, *W_raw; \
  phist_lidx ldV, ldW; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))W,&W_raw,&ldW,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nW_ = std::min(ldW,((nW-1)/cl_size+1)*cl_size); \
  if( nW_+cl_size > ldW ) nW_ = ldW; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,(V!=W),ldV,ldW,*iflag,0, STREAM_TRIAD((nV_+(V!=W)*nW_+nW_)*n*sizeof(_ST_)),4*n*nV*nV);

#endif



//! \def PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(a,V,b,W,iflag)
//! checks performance of mvec_times_sdMat
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(a,V,b,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=_ST_(0)),(b!=_ST_(0)),nV,nW,*iflag,0,0, STREAM_TRIAD(((a!=_ST_(0))*nV+(b!=_ST_(0))*nW+nW)*n*sizeof(_ST_)),2*n*nV*nW);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(a,V,b,W,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag); \
  _ST_ *V_raw, *W_raw; \
  phist_lidx ldV, ldW; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))W,&W_raw,&ldW,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nW_ = std::min(ldW,((nW-1)/cl_size+1)*cl_size); \
  if( nW_+cl_size > ldW ) nW_ = ldW; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,(a!=_ST_(0)),(b!=_ST_(0)),nV,nW,ldV,ldW,*iflag, STREAM_TRIAD(((a!=_ST_(0))*nV_+(b!=_ST_(0))*nW_+nW_)*n*sizeof(_ST_)),2*n*nV*nW);

#endif


//! \def PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag)
//! checks performance of mvec_times_sdMat_inplace
#if !defined(PHIST_PERFCHECK_REALISTIC) || !defined(PHIST_MVECS_ROW_MAJOR)

// ideal model
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&nW,iflag),*iflag); \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,*iflag,0,0,0,0, STREAM_TRIAD((nV+nW)*n*sizeof(_ST_)),n*nV*nW*2);

#else

// realistic model which respects cache line length
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag) \
  int tmp_iflag = *iflag; \
  PHIST_PERFCHECK_MVEC_LEN_T n; \
  int nV, nW; \
  PHIST_CHK_IERR(PHIST_PERFCHECK_MVEC_LENGTH(V,&n,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag); \
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&nW,iflag),*iflag); \
  _ST_ *V_raw; \
  phist_lidx ldV; \
  SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&V_raw,&ldV,iflag); \
  phist_lidx cl_size = phist_cacheline_size<_ST_>(); \
  int nV_ = std::min(ldV,((nV-1)/cl_size+1)*cl_size); \
  if( nV_+cl_size > ldV ) nV_ = ldV; \
  int nW_ = std::min(ldV,((nW-1)/cl_size+1)*cl_size); \
  if( nW_+cl_size > ldV ) nW_ = ldV; \
  *iflag = tmp_iflag; \
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,ldV,*iflag,0,0,0, STREAM_TRIAD((nV_+nW_)*n*sizeof(_ST_)),n*nV*nW*2);

#endif

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
#define PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC_TIMES_SDMAT(V,W,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(a,V,b,W,iflag)
#define PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag)
#define PHIST_PERFCHECK_VERIFY_SPMV(a,A,x,b,y,g,d,iflag)


#endif


#endif /* PHIST_KERNEL_PERFMODELS_HPP */
