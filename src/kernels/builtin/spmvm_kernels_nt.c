/* \file spmvm_kernels_nt.c
 * Fast parallel spMVM kernels with nontemporal stores for different block vector sizes for crsMat_module
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
#include <emmintrin.h>
#include <stdlib.h>

#ifdef PHIST_HIGH_PRECISION_KERNELS
#include "prec_helpers.h"
#else
static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}
#endif


void dspmvm_nt_1_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict shifts, const double *restrict rhsv, const double *restrict halo, double *restrict lhsv)
{
  if( !is_aligned(lhsv,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }
  int nrows2 = nrows/2;

#ifdef PHIST_HIGH_PRECISION_KERNELS
  __m128d shifts_ = _mm_set1_pd(shifts[0]);
  __m128d alpha_ = _mm_set1_pd(alpha);
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows2; i++)
  {
    // apply shifts
    __m128d lhs, lhsC;
    __m128d rhs = _mm_load_pd(rhsv+2*i);
    MM128_2MULTFMA(shifts_,rhs,lhs,lhsC);

    long j, j1;
    for(j = row_ptr[2*i]-1, j1 = row_ptr[2*i+1]-1; j < halo_ptr[2*i]-1 && j1 < halo_ptr[2*i+1]-1; j++, j1++)
    {
      rhs = _mm_set_pd(rhsv[col_idx[j1]-1],rhsv[col_idx[j]-1]);
      __m128d val_ = _mm_set_pd(val[j1],val[j]);
      MM128_4DOTADD(val_,rhs,lhs,lhsC);
    }
    for(;j < halo_ptr[2*i]-1; j++)
    {
      rhs = _mm_set_pd(0.,rhsv[col_idx[j]-1]);
      __m128d val_ = _mm_set_pd(0.,val[j]);
      MM128_4DOTADD(val_,rhs,lhs,lhsC);
    }
    for(;j1 < halo_ptr[2*i+1]-1; j1++)
    {
      rhs = _mm_set_pd(rhsv[col_idx[j1]-1],0.);
      __m128d val_ = _mm_set_pd(val[j1],0.);
      MM128_4DOTADD(val_,rhs,lhs,lhsC);
    }

    for(j = halo_ptr[2*i]-1, j1 = halo_ptr[2*i+1]-1; j < row_ptr[2*i+1]-1 && j1 < row_ptr[2*i+2]-1; j++, j1++)
    {
      rhs = _mm_set_pd(halo[col_idx[j1]-1],halo[col_idx[j]-1]);
      __m128d val_ = _mm_set_pd(val[j1],val[j]);
      MM128_4DOTADD(val_,rhs,lhs,lhsC);
    }
    for(; j < row_ptr[2*i+1]-1; j++)
    {
      rhs = _mm_set_pd(0.,halo[col_idx[j]-1]);
      __m128d val_ = _mm_set_pd(0.,val[j]);
      MM128_4DOTADD(val_,rhs,lhs,lhsC);
    }
    for(; j1 < row_ptr[2*i+2]-1; j1++)
    {
      rhs = _mm_set_pd(halo[col_idx[j1]-1],0.);
      __m128d val_ = _mm_set_pd(val[j1],0.);
      MM128_4DOTADD(val_,rhs,lhs,lhsC);
    }

    _mm_stream_pd(lhsv+2*i,_mm_mul_pd(_mm_add_pd(lhs,lhsC),alpha_));
  }
#else
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows2; i++)
  {
    double lhs1 = shifts[0]*rhsv[2*i];

    for(long j = row_ptr[2*i]-1; j < halo_ptr[2*i]-1; j++)
      lhs1 += val[j]*rhsv[ (col_idx[j]-1) ];
    for(long j = halo_ptr[2*i]-1; j < row_ptr[2*i+1]-1; j++)
      lhs1 += val[j]*halo[ (col_idx[j]-1) ];

    double lhs2 = shifts[0]*rhsv[2*i+1];

    for(long j = row_ptr[2*i+1]-1; j < halo_ptr[2*i+1]-1; j++)
      lhs2 += val[j]*rhsv[ (col_idx[j]-1) ];
    for(long j = halo_ptr[2*i+1]-1; j < row_ptr[2*i+2]-1; j++)
      lhs2 += val[j]*halo[ (col_idx[j]-1) ];

    // multiply with alpha
    lhs1 *= alpha;
    lhs2 *= alpha;

    // non-temporal store of 2 elements
    __m128d lhs_ = _mm_set_pd(lhs2,lhs1);
    _mm_stream_pd(lhsv+2*i, lhs_);
  }
#endif

  // last row
#ifdef PHIST_HIGH_PRECISION_KERNELS
  if( nrows % 2 != 0 )
  {
    double lhs, lhsC;
    DOUBLE_2MULTFMA(shifts[0],rhsv[nrows-1],lhs,lhsC);

    for(long j = row_ptr[nrows-1]-1; j < halo_ptr[nrows-1]-1; j++)
      DOUBLE_4DOTADD(val[j],rhsv[(col_idx[j]-1)],lhs,lhsC);
    for(long j = halo_ptr[nrows-1]-1; j < row_ptr[nrows]-1; j++)
      DOUBLE_4DOTADD(val[j],halo[(col_idx[j]-1)],lhs,lhsC);

    lhsv[nrows-1] = alpha*(lhs+lhsC);
  }
#else
  if( nrows % 2 != 0 )
  {
    lhsv[nrows-1] = shifts[0]*rhsv[nrows-1];
    for(long j = row_ptr[nrows-1]-1; j < halo_ptr[nrows-1]-1; j++)
      lhsv[nrows-1] += val[j]*rhsv[ (col_idx[j]-1) ];
    for(long j = halo_ptr[nrows-1]-1; j < row_ptr[nrows]-1; j++)
      lhsv[nrows-1] += val[j]*halo[ (col_idx[j]-1) ];
    lhsv[nrows-1] *= alpha;
  }
#endif
}


void dspmvm_nt_2_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict shifts, const double *restrict rhsv, const double *restrict halo, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,32) || ldl % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }

  if( !is_aligned(rhsv,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)rhsv);
    exit(1);
  }

  if( !is_aligned(halo,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)halo);
    exit(1);
  }



  __m128d shifts_ = _mm_loadu_pd(shifts);
  __m128d alpha_ = _mm_set1_pd(alpha);

#ifdef PHIST_HIGH_PRECISION_KERNELS
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d rhs_, val_, lhs, lhsC;
    rhs_ = _mm_load_pd(rhsv+2*i);
    MM128_2MULTFMA(shifts_,rhs_,lhs,lhsC);

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      val_ = _mm_load1_pd(val+j);
      rhs_ = _mm_load_pd(rhsv+(col_idx[j]-1)*2);
      MM128_4DOTADD(val_,rhs_,lhs,lhsC);
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      val_ = _mm_load1_pd(val+j);
      rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*2);
      MM128_4DOTADD(val_,rhs_,lhs,lhsC);
    }

    // multiply with alpha
    lhs = _mm_mul_pd(alpha_,_mm_add_pd(lhs,lhsC));
 
    // non-temporal store
    _mm_stream_pd(lhsv+i*ldl, lhs);
  }
#else
#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_ = _mm_load_pd(rhsv+2*i);
    lhs_ = _mm_mul_pd(lhs_, shifts_);

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);

      __m128d rhs_ = _mm_load_pd(rhsv+(col_idx[j]-1)*2);
      rhs_ = _mm_mul_pd(val_,rhs_);
      lhs_ = _mm_add_pd(lhs_,rhs_);
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);

      __m128d rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*2);
      rhs_ = _mm_mul_pd(val_,rhs_);
      lhs_ = _mm_add_pd(lhs_,rhs_);
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    lhs_ = _mm_mul_pd(alpha_,lhs_);
 
    // non-temporal store
    _mm_stream_pd(lhsv+i*ldl, lhs_);
  }
#endif
}

void dspmvm_nt_4_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict shifts, const double *restrict rhsv, const double *restrict halo, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,32) || ldl % 4 != 0 )
  {
    printf("%s: lhsv not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }

  if( !is_aligned(rhsv,32) )
  {
    printf("%s: rhsv not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)rhsv);
    exit(1);
  }

  if( !is_aligned(halo,32) )
  {
    printf("%s: halo not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)halo);
    exit(1);
  }


#ifdef PHIST_HIGH_PRECISION_KERNELS

  __m256d shifts_ = _mm256_loadu_pd(shifts);
  __m256d alpha_ = _mm256_set1_pd(alpha);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m256d rhs_, val_, lhs, lhsC;
    rhs_ = _mm256_load_pd(rhsv+4*i);
    MM256_2MULTFMA(shifts_,rhs_,lhs,lhsC);

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      val_ = _mm256_set1_pd(val[j]);
      rhs_ = _mm256_load_pd(rhsv+(col_idx[j]-1)*4);
      MM256_4DOTADD(val_,rhs_,lhs,lhsC);
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      val_ = _mm256_set1_pd(val[j]);
      rhs_ = _mm256_load_pd(halo+(col_idx[j]-1)*4);
      MM256_4DOTADD(val_,rhs_,lhs,lhsC);
    }

    // multiply with alpha
    lhs = _mm256_mul_pd(alpha_,_mm256_add_pd(lhs,lhsC));
 
    // non-temporal store
    _mm256_stream_pd(lhsv+i*ldl, lhs);
  }

#else

  __m128d shifts_[2];
  shifts_[0] = _mm_loadu_pd(shifts);
  shifts_[1] = _mm_loadu_pd(shifts+2);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[2];
    for(int k = 0; k < 2; k++)
    {
      lhs_[k] = _mm_load_pd(rhsv+4*i+2*k);
      lhs_[k] = _mm_mul_pd(lhs_[k],shifts_[k]);
    }

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 2; k++)
      {
        __m128d rhs_ = _mm_load_pd(rhsv+(col_idx[j]-1)*4+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 2; k++)
      {
        __m128d rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*4+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha and non-temporal store
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 2; k++)
    {
      lhs_[k] = _mm_mul_pd(lhs_[k], alpha_);
      _mm_stream_pd(lhsv+i*ldl+2*k, lhs_[k]);
    }
  }

#endif
}


void dspmvm_nt_8_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                 const double *restrict shifts, const double *restrict rhsv, const double *restrict halo, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) || ldl % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }

  if( !is_aligned(rhsv,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)rhsv);
    exit(1);
  }

  if( !is_aligned(halo,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)halo);
    exit(1);
  }


  __m128d shifts_[4];
  for(int k = 0; k < 4; k++)
    shifts_[k] = _mm_loadu_pd(shifts+2*k);


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[4];
    for(int k = 0; k < 4; k++)
    {
      lhs_[k] = _mm_load_pd(rhsv+8*i+2*k);
      lhs_[k] = _mm_mul_pd(lhs_[k],shifts_[k]);
    }

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 4; k++)
      {
        __m128d rhs_ = _mm_load_pd(rhsv+(col_idx[j]-1)*8+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 4; k++)
      {
        __m128d rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*8+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha and non-temporal store
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 4; k++)
    {
      lhs_[k] = _mm_mul_pd(lhs_[k], alpha_);
      _mm_stream_pd(lhsv+i*ldl+2*k, lhs_[k]);
    }
  }
}


void dspmvm_nt_strided_2_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                         const double *restrict shifts, const double *restrict rhsv, int ldr, const double *restrict halo, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) || ldl % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }

  if( !is_aligned(halo,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)halo);
    exit(1);
  }


  __m128d shifts_ = _mm_loadu_pd(shifts);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_ = _mm_loadu_pd(rhsv+ldr*i);
    lhs_ = _mm_mul_pd(lhs_, shifts_);

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);

      __m128d rhs_ = _mm_loadu_pd(rhsv+(col_idx[j]-1)*ldr);
      rhs_ = _mm_mul_pd(val_,rhs_);
      lhs_ = _mm_add_pd(lhs_,rhs_);
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);

      __m128d rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*2);
      rhs_ = _mm_mul_pd(val_,rhs_);
      lhs_ = _mm_add_pd(lhs_,rhs_);
    }

    // multiply with alpha
    __m128d alpha_ = _mm_set1_pd(alpha);
    lhs_ = _mm_mul_pd(alpha_,lhs_);
 
    // non-temporal store
    _mm_stream_pd(lhsv+i*ldl, lhs_);
  }
}

void dspmvm_nt_strided_4_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                         const double *restrict shifts, const double *restrict rhsv, int ldr, const double *restrict halo, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) || ldl % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }

  if( !is_aligned(halo,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)halo);
    exit(1);
  }


  __m128d shifts_[2];
  shifts_[0] = _mm_loadu_pd(shifts);
  shifts_[1] = _mm_loadu_pd(shifts+2);

#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[2];
    for(int k = 0; k < 2; k++)
    {
      lhs_[k] = _mm_loadu_pd(rhsv+ldr*i+2*k);
      lhs_[k] = _mm_mul_pd(lhs_[k],shifts_[k]);
    }

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 2; k++)
      {
        __m128d rhs_ = _mm_loadu_pd(rhsv+(col_idx[j]-1)*ldr+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 2; k++)
      {
        __m128d rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*4+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha and non-temporal store
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 2; k++)
    {
      lhs_[k] = _mm_mul_pd(lhs_[k], alpha_);
      _mm_stream_pd(lhsv+i*ldl+2*k, lhs_[k]);
    }
  }
}

void dspmvm_nt_strided_8_c(int nrows, double alpha, const long *restrict row_ptr, const long *restrict halo_ptr, const int *restrict col_idx, const double *restrict val,
                         const double *restrict shifts, const double *restrict rhsv, int ldr, const double *restrict halo, double *restrict lhsv, int ldl)
{
  if( !is_aligned(lhsv,16) || ldl % 2 != 0 )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)lhsv);
    exit(1);
  }

  if( !is_aligned(halo,16) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)halo);
    exit(1);
  }


  __m128d shifts_[4];
  for(int k = 0; k < 4; k++)
    shifts_[k] = _mm_loadu_pd(shifts+2*k);


#pragma omp parallel for schedule(static)
  for(int i = 0; i < nrows; i++)
  {
    __m128d lhs_[4];
    for(int k = 0; k < 4; k++)
    {
      lhs_[k] = _mm_loadu_pd(rhsv+ldr*i+2*k);
      lhs_[k] = _mm_mul_pd(lhs_[k],shifts_[k]);
    }

    for(long j = row_ptr[i]-1; j < halo_ptr[i]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 4; k++)
      {
        __m128d rhs_ = _mm_loadu_pd(rhsv+(col_idx[j]-1)*ldr+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    for(long j = halo_ptr[i]-1; j < row_ptr[i+1]-1; j++)
    {
      __m128d val_ = _mm_load1_pd(val+j);
      for(int k = 0; k < 4; k++)
      {
        __m128d rhs_ = _mm_load_pd(halo+(col_idx[j]-1)*8+2*k);
        rhs_ = _mm_mul_pd(val_,rhs_);
        lhs_[k] = _mm_add_pd(lhs_[k],rhs_);
      }
    }

    // multiply with alpha and non-temporal store
    __m128d alpha_ = _mm_set1_pd(alpha);
    for(int k = 0; k < 4; k++)
    {
      lhs_[k] = _mm_mul_pd(lhs_[k], alpha_);
      _mm_stream_pd(lhsv+i*ldl+2*k, lhs_[k]);
    }
  }
}


