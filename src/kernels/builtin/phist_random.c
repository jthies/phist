/*! \file phist_random.c
 * a variant of George Marsaglia's KISS random number generator with jump ahead
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * inspired by https://github.com/cmcqueen/simplerandom and http://web.mst.edu/~vojtat/class_5403/kiss05/rkiss05.f90
 * and http://www.thecodingforums.com/threads/64-bit-kiss-rngs.673657/
 *
*/

#include "phist_random.h"

#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>

#ifdef PHIST_HIGH_PRECISION_KERNELS
#include <immintrin.h>
#endif

//! global state of the random number generator
// (this could be part of a phist_comm, but as we usually only do not use different comms, we can simply define it here)
typedef struct KissState
{
  uint64_t x, y, z, c;
} KissState;

static KissState random_state;


//! initialize the random number generator
void phist_random_init()
{
#pragma omp critical (phist_random)
  {
    random_state.x = 1234567890987654321ULL;
    random_state.y = 362436362436362436ULL;
    random_state.z = 1066149217761810ULL;
    random_state.c = 123456123456123456ULL;
  }
}


// multiply with carry
static inline uint64_t MWC(uint64_t *restrict x, uint64_t *restrict c)
{
  uint64_t t = (*x << 58) + *c;
  *c = *x >> 6;
  *x += t;
  *c += (*x < t);
  return *x;
}

// xorshift
static inline uint64_t XSH(uint64_t *restrict y)
{
  *y ^= (*y<<13);
  *y ^= (*y>>17);
  *y ^= (*y<<43);
  return *y;
}

// congruential
static inline uint64_t CNG(uint64_t *restrict z)
{
  *z = 6906969069ULL * (*z) + 1234567ULL;
  return *z;
}

// KISS transformed to double [-1,1)
static inline double KISSD(uint64_t *restrict x, uint64_t *restrict y, uint64_t *restrict z, uint64_t *restrict c)
{
  return 1.0d/(1ULL<<63) * (MWC(x,c) + XSH(y) + CNG(z)) - 1.0d;
}

//! get a set of random numbers between -1 and +1
void phist_Drandom_number(int n, double*x)
{
#pragma omp critical (phist_random)
  for(int i = 0; i < n; i++)
  {
    x[i] = KISSD(&random_state.x, &random_state.y, &random_state.z, &random_state.c);
  }
}

// skip a list of random numbers (state jump ahead)
// TODO: more sophisticated O(log n) code
static inline void KISS_SKIP(uint64_t n, uint64_t *restrict x, uint64_t *restrict y, uint64_t *restrict z, uint64_t *restrict c)
{
  for(int i = 0; i < n; i++)
  {
    KISSD(x,y,z,c);
  }
}

// check alignement
static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count)
{
  return (uintptr_t)pointer % byte_count == 0;
}


// generate random numbers in parallel for a dense block
// TODO: 4-way unrolling of KISS state using leapfrogging for performance
void drandom_1(int nrows, double *restrict v, int64_t pre_skip, int64_t post_skip)
{
  if( !is_aligned(v,32) )
  {
    printf("%s: not aligned %lx\n", __FUNCTION__, (uintptr_t)(void*)v);
    exit(1);
  }

#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

#pragma omp critical (phist_random)
  {
    // calculate array offsets per thread (handcoded parallel for schedule(static))
    int threadOffset4[nt+1];
    int nrows4 = nrows/4;
    int rows4PerThread = nrows4/nt;
    int rows4ModThread = nrows4%nt;
    threadOffset4[0] = 0;
    for(int i = 1; i <= rows4ModThread; i++)
      threadOffset4[i] = threadOffset4[i-1] + rows4PerThread+1;
    for(int i = rows4ModThread+1; i <= nt; i++)
      threadOffset4[i] = threadOffset4[i-1] + rows4PerThread;

    uint64_t x = random_state.x;
    uint64_t y = random_state.y;
    uint64_t z = random_state.z;
    uint64_t c = random_state.c;
#pragma omp parallel firstprivate(x,y,z,c)
    {
#ifdef PHIST_HAVE_OPENMP
      int it = omp_get_thread_num();
#else
      int it = 1;
#endif

      // jump to first index of this thread
      KISS_SKIP(pre_skip+4*threadOffset4[it], &x,&y,&z,&c);

      // loop
      for(int i = threadOffset4[it]; i < threadOffset4[it+1]; i++)
      {
        double yi1 = KISSD(&x,&y,&z,&c);
        double yi2 = KISSD(&x,&y,&z,&c);
        double yi3 = KISSD(&x,&y,&z,&c);
        double yi4 = KISSD(&x,&y,&z,&c);
#ifdef PHIST_HIGH_PRECISION_KERNELS
        // use AVX non-temporal stores
        __m256d yi = _mm256_set_pd(yi4,yi3,yi2,yi1);
        _mm256_stream_pd(v+4*i, yi);
#else
        // use SSE non-temporal stores
        __m128d yi21 = _mm_set_pd(yi2,yi1);
        _mm_stream_pd(v+4*i, yi21);
        __m128d yi43 = _mm_set_pd(yi4,yi3);
        _mm_stream_pd(v+4*i+2, yi43);
#endif
      }

      // the last thread cares about the remainder and setting the state
      if( it == nt-1 )
      {
        for(int i = 4*threadOffset4[it+1]; i < nrows; i++)
          v[i] = KISSD(&x,&y,&z,&c);

        // post skip (to be in sync with other MPI procs
        KISS_SKIP(post_skip, &x,&y,&z,&c);
        random_state.x = x;
        random_state.y = y;
        random_state.z = z;
        random_state.c = c;
      }
    }
  }
}


// generate random numbers in parallel for strided block
// TODO: 4-way unrolling of KISS state using leapfrogging for performance
void drandom_general(int nvec, int nrows, double *restrict v, int ldv, int64_t pre_skip, int64_t post_skip)
{
#ifdef PHIST_HAVE_OPENMP
  int nt = omp_get_max_threads();
#else
  int nt = 1;
#endif

#pragma omp critical (phist_random)
  {
    // calculate array offsets per thread (handcoded parallel for schedule(static))
    int threadOffset[nt+1];
    int rowsPerThread = nrows/nt;
    int rowsModThread = nrows%nt;
    threadOffset[0] = 0;
    for(int i = 1; i <= rowsModThread; i++)
      threadOffset[i] = threadOffset[i-1] + rowsPerThread+1;
    for(int i = rowsModThread+1; i <= nt; i++)
      threadOffset[i] = threadOffset[i-1] + rowsPerThread;

    uint64_t x = random_state.x;
    uint64_t y = random_state.y;
    uint64_t z = random_state.z;
    uint64_t c = random_state.c;
#pragma omp parallel firstprivate(x,y,z,c)
    {
#ifdef PHIST_HAVE_OPENMP
      int it = omp_get_thread_num();
#else
      int it = 1;
#endif

      // jump to first index of this thread
      KISS_SKIP(pre_skip+threadOffset[it], &x,&y,&z,&c);

      // loop
      for(int i = threadOffset[it]; i < threadOffset[it+1]; i++)
        for(int j = 0; j < nvec; j++)
          v[i*ldv+j] = KISSD(&x,&y,&z,&c);

      // the last thread cares about the remainder and setting the state
      if( it == nt-1 )
      {
        // post skip (to be in sync with other MPI procs
        KISS_SKIP(post_skip, &x,&y,&z,&c);
        random_state.x = x;
        random_state.y = y;
        random_state.z = z;
        random_state.c = c;
      }
    }
  }
}

