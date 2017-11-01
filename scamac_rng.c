#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "scamac_rng.h"

scamac_rng_st * scamac_rng_alloc(uint64_t seed) {
  scamac_rng_st * rg = malloc(sizeof *rg);
  scamac_rng_reset(seed, rg);
  return rg;
}

void scamac_rng_free(scamac_rng_st * rg) {
  if (rg) {
    free(rg);
  }
}

void scamac_rng_reset(uint64_t seed, scamac_rng_st *rg) {
  int i;

  rg->seed=seed;
  for (i=0; i<16; i++) {
    rg->sinit[i]=seed;
    rg->s[i]=rg->sinit[i];
  }

  //warm up
  rg->p=0;
  for (i=0; i<1000; i++) {
    scamac_rng_next(rg);
  }
  rg->idx=0;
}

// xorshift random number generator according to Sebastiano Vigna, arxiv:1402.6246
// [ ACM Trans. Math. Software 42, 4 (2016) ]
uint64_t scamac_rng_next(scamac_rng_st *rg) {
  uint64_t s0 = rg->s[rg->p];
  rg->p = (rg->p + 1) & 15;
  uint64_t s1 = rg->s[rg->p];
  s1 ^= s1 << 31; // a
  rg->s[rg->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c

  rg->prn = rg->s[rg->p] * UINT64_C(1181783497276652981);
  (rg->idx)++;

  return rg->prn;
}

uint64_t scamac_rng_get(scamac_rng_st * rg, int64_t i) {
  if (rg->idx > i) {
    scamac_rng_reset(rg->seed,rg);
  }
  while (rg->idx < i) {
    scamac_rng_next(rg);
  }
  return rg->prn;
}

// uniform in [a,b)
double scamac_rng_get_double(scamac_rng_st *rg, double a, double b, int64_t i) {
  uint64_t ui = scamac_rng_get(rg,i);
  unsigned long int u1 = ui & 0xFFFFFFFF;
  unsigned long int u2 = (ui >> 32) & 0xFFFFFFFF;
  double v = ldexp((double) u1,-64) + ldexp((double) u2,-32);
  return a + (b-a)*v;
}

uint64_t scamac_generate_seed(int short_seed) {
  int ti = (int)time(NULL);  // do whatever you want here
  return (uint64_t) short_seed + ( ((uint64_t) ti) << 32 );
}
