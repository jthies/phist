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
#include <string.h>

#ifdef PHIST_HIGH_PRECISION_KERNELS
#include <immintrin.h>
#endif


// ************************************************************************
//                            GLOBAL STATE
// ************************************************************************

//! global state of the random number generator
// (this could be part of a phist_comm, but as we usually only do not use different comms, we can simply define it here)
typedef struct KissState
{
  uint64_t x, y, z, c;
} KissState;

//! the global state, should be identical on all MPI processes (intrinsically!)
static KissState random_state;

//! flag to check if skip ahead lookup tables are initialized
static _Bool lookup_tables_initialized = 0;



// ************************************************************************
//                            PUBLIC FUNCTIONS
// ************************************************************************

// declarations, description see implementation below
static void init_lookup_tables();
static inline double KISSD(uint64_t *restrict x, uint64_t *restrict y, uint64_t *restrict z, uint64_t *restrict c);
static inline void KISS_SKIP(uint64_t n, uint64_t *restrict x, uint64_t *restrict y, uint64_t *restrict z, uint64_t *restrict c);
static inline _Bool is_aligned(const void *restrict pointer, size_t byte_count) {return (uintptr_t)pointer % byte_count == 0; }


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

  if( !lookup_tables_initialized )
    init_lookup_tables();

  lookup_tables_initialized = 1;
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
      KISS_SKIP(pre_skip+nvec*threadOffset[it], &x,&y,&z,&c);

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



// ************************************************************************
//                            SKIP AHEAD LOOKUP TABLES
// ************************************************************************

// tables for XSH and XSH_SKIP
static const int_fast8_t XSH_shift[3] = {13,-17,43};
static uint64_t XSH_matrixpow[64][64];

// tables for CNG and CNG_SKIP
static const uint64_t CNG_r = 6906969069ULL;
static const uint64_t CNG_c = 1234567ULL;
static uint64_t CNG_powr[128];

// tables for MWC and MWC_SKIP
// MWC with factor 2^58+1 and modulus 2^64
// is equivalent to LCG with modulus m = (2^58+1)*2^64-1
static uint64_t MWC_lcg_m[2]   = {(uint64_t)-1, (uint64_t)1 << 58};
// and the factor is the multiplicative inverse of 2^64 modulus m
static uint64_t MWC_lcg_a[2] = {288230376151711745ULL, 0ULL};
static uint64_t MWC_lcg_powa[64][2];



// ************************************************************************
//                            SMALL HELPER FUNCTIONS
// ************************************************************************

// XSH_SKIP helper function to generate a shift matrix
static void bitmatrix_I_shift(uint64_t M[64], int_fast8_t shift)
{
  uint64_t eye = 1ULL;
  uint64_t val = shift >= 0 ? 1ULL<<shift : 0ULL;
  for(int i = 0; i < 64; i++)
  {
    M[i] = eye ^ val;
    eye <<= 1ULL;
    if( shift < 0 )
    {
      shift++;
      if( shift == 0 )
        val = 1ULL;
    }
    else
    {
      val <<= 1ULL;
    }
  }
}

// XSH_SKIP helper function to apply a bitmatrix to a vector
static uint64_t bitmatrix_apply_vec(const uint64_t M[64], uint64_t v)
{
  uint64_t r = 0ULL;
  for(int i = 0; i < 64; i++)
  {
    if( v & 1ULL )
      r ^= M[i];
    v >>= 1;
  }
  return r;
}

// XSH_SKIP helper function to apply a bitmatrix to another matrix
// result is stored in M, e.g. M *= N
static void bitmatrix_apply_mat(uint64_t M[64], const uint64_t N[64])
{
  uint64_t tmp[64];
  memcpy(tmp, M, 64*sizeof(uint64_t));
  for(int i = 0; i < 64; i++)
    M[i] = bitmatrix_apply_vec(tmp, N[i]);
}


// MWC_SKIP helper function for shifting unsigned 128bit integer (by at most 64)
static inline void lshift(uint64_t x[2], const uint_fast8_t s)
{
  uint64_t t = x[0] >> (64-s);
  x[1] <<= s;
  x[1] ^= t;
  x[0] <<= s;
}

// MWC_SKIP helper function for shifting unsigned 128bit integer (by at most 64)
static inline void rshift(uint64_t x[2], const uint_fast8_t s)
{
  uint64_t t = x[1] << (64-s);
  x[0] >>= s;
  x[0] ^= t;
  x[1] >>= s;
}

// MWC_SKIP helper function for addition of two unsigned 128 bit integers
static inline void add(uint64_t a[2], const uint64_t b[2])
{
  a[0] += b[0];
  a[1] += b[1] + (a[0]<b[0]);
}

// MWC_SKIP helper function for subtraction of two unsigned 128 bit integers
static inline void sub(uint64_t a[2], const uint64_t b[2])
{
  a[1] -= b[1] + (a[0]<b[0]);
  a[0] -= b[0];
}

// MWC_SKIP helper function for the comparison of two unsigned 128 bit integers
static inline _Bool ge(const uint64_t a[2], const uint64_t b[2])
{
  return a[1] > b[1] || (a[1] == b[1] && a[0] >= b[0]);
}

// MWC_SKIP helper function to calculate usigned integer 128bit modulus MWC_lcg_m
static void mod_MWC(uint64_t a[2])
{
  for(int i = 5; i >= 0; i--)
  {
    uint64_t modp2[2] = {MWC_lcg_m[0], MWC_lcg_m[1]};
    lshift(modp2,i);
    if( ge(a,modp2) )
      sub(a,modp2);
  }
}

// MWC_SKIP helper function to multiply two 128bit numbers modulo MWC_lcg_m
static void mul_mod_MWC(uint64_t a[2], const uint64_t b[2])
{
  uint64_t b_[2] = {b[0],b[1]};
  mod_MWC(b_);

  uint64_t a_[2] = {a[0],a[1]};
  mod_MWC(a_);

  a[1] = 0ULL; a[0] = 0ULL;
  for(int i = 0; i < 128; i++)
  {
    if( a_[0] & 1ULL )
    {
      add(a,b_);
      if( ge(a,MWC_lcg_m) )
        sub(a,MWC_lcg_m);
    }

    rshift(a_,1);
    if( a_[1] == 0 && a_[0] == 0 )
      break;

    lshift(b_,1);
    if( ge(b_,MWC_lcg_m) )
      sub(b_,MWC_lcg_m);
  }
}



// ************************************************************************
//                            RANDOM NUMBER GENERATORS
// ************************************************************************

// multiply with carry
static inline uint64_t MWC(uint64_t *restrict x, uint64_t *restrict c)
{
  // this actually does:
  // newX = ((2^58+1)*x + c) mod 2^64
  // newC = ((2^58+1)*x + c) / 2^64    // rounded down
  // e.g. this is equivalent to the following calculation where y/newY can be represented as uin128_t
  // (but be aware of overflows!)
  // y = 2^64*c + x   (mod 2^128)
  // newY = MWC_lcg_a * y mod MWC_lcg_m
  // where MWC_lcg_m = (2^58+1)*2^64-1
  // and   MWC_lcg_a = (2^64)^(-1) mod m
  // newX = (uint64_t) newY;          // lower part
  // newC = (uint64_t) (newY >> 64);  // upper part
  uint64_t t = (*x << 58) + *c;
  *c = *x >> 6;
  *x += t;
  *c += (*x < t);
  /*
   * equivalent to:
    uint64_t t[2] = {*x,*c};
    mul_mod_MWC(t,MWC_lcg_a);
    *x = t[0];
    *c = t[1];
  */
  return *x;
}

// multiply with carry skip
static inline void MWC_SKIP(uint64_t n, uint64_t *restrict x, uint64_t *restrict c)
{
  // exploit
  // y_k = (a^k y) mod m = (a^k mod m) * y mod m

  // calculate a^k mod m:
  uint64_t a[2] = {1ULL,0ULL};
  int i = 0;
  for(; i < 64; i++)
  {
    if( n & 1ULL )
    {
      a[0] = MWC_lcg_powa[i][0];
      a[1] = MWC_lcg_powa[i][1];
      i++;
      n >>= 1;
      break;
    }
    n >>= 1;
  }
  for(; i < 64; i++)
  {
    if( n & 1ULL )
      mul_mod_MWC(a, MWC_lcg_powa[i]);
    n >>= 1;
    if( n == 0 )
      break;
  }

  // calculate a*(c,x)
  uint64_t t[2] = {*x,*c};
  mul_mod_MWC(t, a);
  *x = t[0];
  *c = t[1];
}


// xorshift
static inline uint64_t XSH(uint64_t *restrict y)
{
  // XSH_shift
  *y ^= (*y<<13);
  *y ^= (*y>>17);
  *y ^= (*y<<43);
  return *y;
}

// xorshift skip
static inline void XSH_SKIP(uint64_t n, uint64_t *restrict y)
{
  // idea: represent shift as a 64x64bit matrix and calculate matrix power
  // calculate XSH_matrix ^ n
  // quite costly, so simply apply XSH for small values of n...
  if( n < 1024 )
  {
    for(int i = 0; i < n; i++)
      XSH(y);
    return;
  }
  uint64_t M[64];
  int i = 0;
  for(; i < 64; i++)
  {
    if( n & 1ULL )
    {
      memcpy(M, XSH_matrixpow[i], 64*sizeof(int64_t));
      i++;
      n >>= 1;
      break;
    }
    n >>= 1;
  }
  for(; i < 64; i++)
  {
    if( n & 1ULL )
      bitmatrix_apply_mat(M, XSH_matrixpow[i]);
    n >>= 1;
    if( n == 0 )
      break;
  }

  // apply matrix
  *y = bitmatrix_apply_vec(M, *y);
}


// congruential
static inline uint64_t CNG(uint64_t *restrict z)
{
  *z = CNG_r * (*z) + CNG_c;
  return *z;
}

// congruential skip
// based on: z * (r^n) + c * (1 + r + r^2 + ... + r^(n-1))  (mod 2^64)
static inline void CNG_SKIP(uint64_t n, uint64_t *restrict z)
{
  // calculate: r^n (mod 2^64)
  uint64_t r = 1ULL;
  uint64_t n_ = n;
  for(int i = 0; i < 64; i++)
  {
    if( n_ & 1ULL )
      r *= CNG_powr[i];
    n_ >>= 1;
    if( n_ == 0 )
      break;
  }

  // calculate geometric series: 1+r^2+...+r^(n-1) (mod 2^64)
  uint64_t c = 0ULL;
  uint64_t tmp = CNG_c;
  n_ = n;
  for(int i = 0; i < 64; i++)
  {
    if( n_ & 1ULL )
    {
      // calculate CNG_powr[i] ^(n-1)
      uint64_t powr = 1ULL;
      uint64_t n__ = n_-1;
      for(int i_ = i; i_ < i+64; i_++)
      {
        if( n__ & 1ULL )
          powr *= CNG_powr[i_];
        n__ >>= 1;
        if( n__ == 0 )
          break;
      }
      c += tmp * powr;
    }
    tmp *= (1 + CNG_powr[i]);
    n_ >>= 1;
    if( n_ == 1 )
      break;
  }
  c += tmp;

  *z = r * *z + c;
}


// KISS transformed to double [-1,1)
static inline double KISSD(uint64_t *restrict x, uint64_t *restrict y, uint64_t *restrict z, uint64_t *restrict c)
{
  return 1.0/(1ULL<<63) * (MWC(x,c) + XSH(y) + CNG(z)) - 1.0;
}

// skip a list of random numbers (state jump ahead)
static inline void KISS_SKIP(uint64_t n, uint64_t *restrict x, uint64_t *restrict y, uint64_t *restrict z, uint64_t *restrict c)
{
  MWC_SKIP(n,x,c);
  XSH_SKIP(n,y);
  CNG_SKIP(n,z);
}


// initialize skip ahead lookup tables
static void init_lookup_tables()
{
  // generate bit matrix for XSH
  // y ^= (y<<13);
  // e.g. y <- (I+shift(13)) * y
  // y ^= (y>>17);
  // e.g. y <- (I+shift(-17)) * y
  // y ^= (y<<43);
  // e.g. y <- (I+shift(43)) * y
  // => (I+shift(43)) * (I+shift(-17)) * (I+shift(13))
  bitmatrix_I_shift(XSH_matrixpow[0], XSH_shift[2]);
  uint64_t tmp[64];
  bitmatrix_I_shift(tmp, XSH_shift[1]);
  bitmatrix_apply_mat(XSH_matrixpow[0], tmp);
  bitmatrix_I_shift(tmp, XSH_shift[0]);
  bitmatrix_apply_mat(XSH_matrixpow[0], tmp);
  // now calculate matrix powers of 2
  for(int i = 1; i < 64; i++)
  {
    memcpy(XSH_matrixpow[i], XSH_matrixpow[i-1], 64*sizeof(int64_t));
    bitmatrix_apply_mat(XSH_matrixpow[i], XSH_matrixpow[i]);
  }


  // generate lookup table for CNG
  CNG_powr[0] = CNG_r;
  for(int i = 1; i < 128; i++)
    CNG_powr[i] = CNG_powr[i-1]*CNG_powr[i-1];


#ifdef TESTING
  // veryify MWC_lcg_a*2^64 = 1 mod MWC_lcg_m
  uint64_t t[2] = {MWC_lcg_a[0], MWC_lcg_a[1]};
  uint64_t b[2] = {1ULL,0ULL};
  mul_mod_MWC(t,b);
  if( t[0] != 1ULL && t[1] != 0ULL )
  {
    printf("phist_random selftest: lcg_a * 2^64 mod lcg_m: %llu,%llu\n", t[1], t[0]);
    exit(1);
  }
#endif
  // generate lookup table for MWC
  MWC_lcg_powa[0][0] = MWC_lcg_a[0];
  MWC_lcg_powa[0][1] = MWC_lcg_a[1];
  for(int i = 1; i < 64; i++)
  {
    MWC_lcg_powa[i][0] = MWC_lcg_powa[i-1][0];
    MWC_lcg_powa[i][1] = MWC_lcg_powa[i-1][1];
    mul_mod_MWC(MWC_lcg_powa[i], MWC_lcg_powa[i-1]);
  }


#ifdef TESTING
  // check result
  uint64_t x = random_state.x;
  uint64_t y = random_state.y;
  uint64_t z = random_state.z;
  uint64_t c = random_state.c;
  for(int i = 0; i < 100000000; i++)
    KISSD(&x,&y,&z,&c);
  uint64_t sum = x+y+z;
  if( sum != 1666297717051644203ULL )
  {
    printf("phist_random selftest: result after 10^8 steps: %llu (should be %llu)\n", x+y+z, 1666297717051644203ULL);
    exit(1);
  }
  // check skipping
  x = random_state.x;
  y = random_state.y;
  z = random_state.z;
  c = random_state.c;
  KISS_SKIP(100000000, &x,&y,&z,&c);
  sum = x+y+z;
  if( sum != 1666297717051644203ULL )
  {
    printf("phist_random selftest: result after skipping 10^8 steps: %llu (should be %llu)\n", x+y+z, 1666297717051644203ULL);
    exit(1);
  }
#endif
}

