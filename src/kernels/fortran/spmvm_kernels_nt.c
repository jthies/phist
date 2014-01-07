// hopefully fast spMVM kernels with nontemporary stores
#include <stdint.h>
#include <stdio.h>
#include <emmintrin.h>

static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}

// provide possibility to check alignment in fortran
void mem_is_aligned16_(const void*restrict pointer, int* ret)
{
  if( is_aligned(pointer,16) )
    *ret = 0;
  else
    *ret = 1;
}

void dspmvm_nt_1_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict rhsv, double *restrict lhsv)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows/2; i++)
  {
    double lhs1 = 0.;

    for(long j = row_ptr[2*i]-1; j < row_ptr[2*i+1]-1; j++)
      lhs1 += val[j]*rhsv[ (col_idx[j]-1) ];

    double lhs2 = 0.;

    for(long j = row_ptr[2*i+1]-1; j < row_ptr[2*i+2]-1; j++)
      lhs2 += val[j]*rhsv[ (col_idx[j]-1) ];

    // multiply with alpha
    lhs1 *= alpha;
    lhs2 *= alpha;

    // non-temporal store of 2 elements
    __m128d lhs_ = _mm_set_pd(lhs2,lhs1);
    _mm_stream_pd(&lhsv[2*i], lhs_);
  }

  // last row
  if( nrows % 2 != 0 )
  {
    lhsv[nrows-1] = 0.;
    for(long j = row_ptr[nrows-1]-1; j < row_ptr[nrows]-1; j++)
      lhsv[nrows-1] += val[j]*rhsv[ (col_idx[j]-1) ];
    lhsv[nrows-1] *= alpha;
  }
}


void dspmvm_nt_2_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict rhsv, double *restrict lhsv, int ldv)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_ = _mm_set1_pd(0.);

    for(long j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_set1_pd(val[j]);

      const double *rhsp = rhsv + (col_idx[j]-1)*2;
      __m128d rhs_ = _mm_set_pd(rhsp[0],rhsp[1]);
      rhs_ = _mm_mul_pd(val_,rhs_);
      lhs_ = _mm_add_pd(lhs_,rhs_);
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    lhs_ = _mm_mul_pd(alpha_,lhs_);
 
    // non-temporal store
    _mm_stream_pd(&lhsv[i*ldv], lhs_);
  }
}

void dspmvm_nt_4_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict rhsv, double *restrict lhsv, int ldv)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[2];
    for(int k = 0; k < 2; k++)
      lhs_[k] = _mm_set1_pd(0.);

    for(long j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_set1_pd(val[j]);
      for(int k = 0; k < 2; k++)
      {
        const double *rhsp = rhsv + (col_idx[j]-1)*4 + 2*k;
        __m128d rhs_ = _mm_set_pd(rhsp[1],rhsp[0]);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 2; k++)
      lhs_[k] = _mm_mul_pd(alpha_,lhs_[k]);

    // non-temporal store
    for(int k = 0; k < 2; k++)
      _mm_stream_pd(&lhsv[i*ldv+2*k], lhs_[k]);
  }
}

void dspmvm_nt_8_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict rhsv, double *restrict lhsv, int ldv)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[4];
    for(int k = 0; k < 4; k++)
      lhs_[k] = _mm_set1_pd(0.);

    for(long j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_set1_pd(val[j]);
      for(int k = 0; k < 4; k++)
      {
        const double *rhsp = rhsv + (col_idx[j]-1)*8 + 2*k;
        __m128d rhs_ = _mm_set_pd(rhsp[1],rhsp[0]);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 4; k++)
      lhs_[k] = _mm_mul_pd(alpha_,lhs_[k]);

    // non-temporal store
    for(int k = 0; k < 4; k++)
      _mm_stream_pd(&lhsv[i*ldv+2*k], lhs_[k]);
  }
}


void dspmvm_nt_strided_2_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                         const double *restrict rhsv, int ldr, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_ = _mm_set1_pd(0.);

    for(long j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_set1_pd(val[j]);

      const double *rhsp = rhsv + (col_idx[j]-1)*ldr;
      __m128d rhs_ = _mm_set_pd(rhsp[0],rhsp[1]);
      rhs_ = _mm_mul_pd(val_,rhs_);
      lhs_ = _mm_add_pd(lhs_,rhs_);
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    lhs_ = _mm_mul_pd(alpha_,lhs_);
 
    // non-temporal store
    _mm_stream_pd(&lhsv[i*ldl], lhs_);
  }
}

void dspmvm_nt_strided_4_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                         const double *restrict rhsv, int ldr, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[2];
    for(int k = 0; k < 2; k++)
      lhs_[k] = _mm_set1_pd(0.);

    for(long j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_set1_pd(val[j]);
      for(int k = 0; k < 2; k++)
      {
        const double *rhsp = rhsv + (col_idx[j]-1)*ldr + 2*k;
        __m128d rhs_ = _mm_set_pd(rhsp[1],rhsp[0]);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 2; k++)
      lhs_[k] = _mm_mul_pd(alpha_,lhs_[k]);

    // non-temporal store
    for(int k = 0; k < 2; k++)
      _mm_stream_pd(&lhsv[i*ldl+2*k], lhs_[k]);
  }
}

void dspmvm_nt_strided_8_c(int nrows, double alpha, const long *restrict row_ptr, const int *restrict col_idx, const double *restrict val,
                         const double *restrict rhsv, int ldr, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("not aligned %lx\n", (uintptr_t)(void*)lhsv);
    exit(1);
  }

#pragma omp parallel for
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[4];
    for(int k = 0; k < 4; k++)
      lhs_[k] = _mm_set1_pd(0.);

    for(long j = row_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_set1_pd(val[j]);
      for(int k = 0; k < 4; k++)
      {
        const double *rhsp = rhsv + (col_idx[j]-1)*ldr + 2*k;
        __m128d rhs_ = _mm_set_pd(rhsp[1],rhsp[0]);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 4; k++)
      lhs_[k] = _mm_mul_pd(alpha_,lhs_[k]);

    // non-temporal store
    for(int k = 0; k < 4; k++)
      _mm_stream_pd(&lhsv[i*ldl+2*k], lhs_[k]);
  }
}


