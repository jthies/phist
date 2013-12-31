// hopefully fast axpby kernels with nontemporary stores
#include <stdint.h>
#include <stdio.h>
#include <emmintrin.h>

static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}


void daxpy_nt_2(int nrows, const double *restrict alpha, const double *restrict x, double *restrict y)
{
  if( !is_aligned(y,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    // get x
    const double *xp = x + 2*i;
    __m128d x_ = _mm_set_pd(xp[1],xp[0]);
    // multiply with alpha
    __m128d alpha_ = _mm_set_pd(alpha[1],alpha[0]);
    __m128d y_ = _mm_mul_pd(x_,alpha_);
    // non-temporal store
    _mm_stream_pd(&y[4*i], y_);
  }
}


void daxpy_nt_4(int nrows, const double *restrict alpha, const double *restrict x, double *restrict y)
{
  if( !is_aligned(y,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    for(int k = 0; k < 2; k++)
    {
      // get x
      const double *xp = x + 4*i + 2*k;
      __m128d x_ = _mm_set_pd(xp[1],xp[0]);
      // multiply with alpha
      __m128d alpha_ = _mm_set_pd(alpha[2*k+1],alpha[2*k]);
      __m128d y_ = _mm_mul_pd(x_,alpha_);
      // non-temporal store
      _mm_stream_pd(&y[4*i+2*k], y_);
    }
  }
}


void daxpy_nt_8(int nrows, const double *restrict alpha, const double *restrict x, double *restrict y)
{
  if( !is_aligned(y,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)y);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    for(int k = 0; k < 4; k++)
    {
      // get x
      const double *xp = x + 8*i + 4*k;
      __m128d x_ = _mm_set_pd(xp[1],xp[0]);
      // multiply with alpha
      __m128d alpha_ = _mm_set_pd(alpha[4*k+1],alpha[4*k]);
      __m128d y_ = _mm_mul_pd(x_,alpha_);
      // non-temporal store
      _mm_stream_pd(&y[8*i+4*k], y_);
    }
  }
}


