/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file axpy_kernels_nt.c
 * Fast parallel BLAS-axpy like functions with nontemporal stores (NT) for different blocksizes for mvec_module
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include <stdint.h>
#include <stdio.h>
#ifdef PHIST_HAVE_SSE
#include <emmintrin.h>
#endif
#include <stdlib.h>

static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}


void daxpy_nt_2_c(int nrows, const double *restrict alpha, const double *restrict x, double *restrict y)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( !is_aligned(y,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

  if( !is_aligned(x,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
  }


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    // get x
    __m128d x_ = _mm_load_pd(x+2*i);
    // multiply with alpha
    __m128d alpha_ = _mm_set_pd(alpha[1],alpha[0]);
    __m128d y_ = _mm_mul_pd(x_,alpha_);
    // non-temporal store
    _mm_stream_pd(y+2*i, y_);
  }
#endif
}


void daxpy_nt_4_c(int nrows, const double *restrict alpha, const double *restrict x, double *restrict y)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( !is_aligned(y,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

  if( !is_aligned(x,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
  }


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    for(int k = 0; k < 2; k++)
    {
      // get x
      __m128d x_ = _mm_load_pd(x+4*i+2*k);
      // multiply with alpha
      __m128d alpha_ = _mm_set_pd(alpha[2*k+1],alpha[2*k]);
      __m128d y_ = _mm_mul_pd(x_,alpha_);
      // non-temporal store
      _mm_stream_pd(y+4*i+2*k, y_);
    }
  }
#endif
}


void daxpy_nt_8_c(int nrows, const double *restrict alpha, const double *restrict x, double *restrict y)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( !is_aligned(y,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

  if( !is_aligned(x,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
  }


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    for(int k = 0; k < 4; k++)
    {
      // get x
      __m128d x_ = _mm_load_pd(x+8*i+2*k);
      // multiply with alpha
      __m128d alpha_ = _mm_set_pd(alpha[2*k+1],alpha[2*k]);
      __m128d y_ = _mm_mul_pd(x_,alpha_);
      // non-temporal store
      _mm_stream_pd(y+8*i+2*k, y_);
    }
  }
#endif
}


void daxpy_nt_strided_2_c(int nrows, const double *restrict alpha, const double *restrict x, int ldx, double *restrict y, int ldy)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( !is_aligned(y,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    // get x
    __m128d x_ = _mm_loadu_pd(x+ldx*i);
    // multiply with alpha
    __m128d alpha_ = _mm_set_pd(alpha[1],alpha[0]);
    __m128d y_ = _mm_mul_pd(x_,alpha_);
    // non-temporal store
    _mm_stream_pd(y+ldy*i, y_);
  }
#endif
}


void daxpy_nt_strided_4_c(int nrows, const double *restrict alpha, const double *restrict x, int ldx, double *restrict y, int ldy)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( !is_aligned(y,16) || ldy % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    for(int k = 0; k < 2; k++)
    {
      // get x
      __m128d x_ = _mm_loadu_pd(x+ldx*i+2*k);
      // multiply with alpha
      __m128d alpha_ = _mm_set_pd(alpha[2*k+1],alpha[2*k]);
      __m128d y_ = _mm_mul_pd(x_,alpha_);
      // non-temporal store
      _mm_stream_pd(y+ldy*i+2*k, y_);
    }
  }
#endif
}


void daxpy_nt_strided_8_c(int nrows, const double *restrict alpha, const double *restrict x, int ldx, double *restrict y, int ldy)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( !is_aligned(y,16) || ldy % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    for(int k = 0; k < 4; k++)
    {
      // get x
      __m128d x_ = _mm_loadu_pd(x+ldx*i+2*k);
      // multiply with alpha
      __m128d alpha_ = _mm_set_pd(alpha[2*k+1],alpha[2*k]);
      __m128d y_ = _mm_mul_pd(x_,alpha_);
      // non-temporal store
      _mm_stream_pd(y+ldy*i+2*k, y_);
    }
  }
#endif
}


void dcopy_general_nt_c(int nrows, int nvec, const double *restrict x, int ldx, double *restrict y, int ldy)
{
#ifndef PHIST_HAVE_SSE
  printf("%s: must not be called on platforms without SSE.", __FUNCTION__);
  exit(1);
#else
  if( nvec % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)x);
    exit(1);
  }

  if( !is_aligned(y,16) || ldy % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    for(int j = 0; j < nvec/2; j++)
    {
      __m128d tmp = _mm_loadu_pd(x+i*ldx+2*j);
      // non-temporal store
      _mm_stream_pd(y+i*ldy+2*j, tmp);
    }
  }
#endif
}

