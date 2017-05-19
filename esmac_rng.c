#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "esmac_rng.h"

esmac_rng_t * esmac_rng_alloc(uint64_t seed) {
  esmac_rng_t * rg = malloc(sizeof *rg);
  esmac_rng_reset(seed, rg);
  return rg;
}

void esmac_rng_free(esmac_rng_t * rg) {
  if (rg) {free(rg);}
}

void esmac_rng_reset(uint64_t seed, esmac_rng_t *rg) {
  int i;
  
  rg->seed=seed;
  for (i=0;i<16;i++) {
    rg->sinit[i]=seed;
    rg->s[i]=rg->sinit[i];
  }
      
  //warm up
  rg->p=0;
  for (i=0;i<1000;i++) {
    esmac_rng_next(rg);
  }
  rg->idx=0;
}

// xorshift random number generator according to Sebastiano Vigna, arxiv:1402.6246
// [ ACM Trans. Math. Software 42, 4 (2016) ]
uint64_t esmac_rng_next(esmac_rng_t *rg) {
  uint64_t s0 = rg->s[rg->p];
  rg->p = (rg->p + 1) & 15;
  uint64_t s1 = rg->s[rg->p];
  s1 ^= s1 << 31; // a
  rg->s[rg->p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
  
  rg->prn = rg->s[rg->p] * UINT64_C(1181783497276652981); 
  (rg->idx)++;
  
  return rg->prn;
}

uint64_t esmac_rng_get(esmac_rng_t * rg, int64_t i) {
  if (rg->idx > i) {
    esmac_rng_reset(rg->seed,rg);
  }
  while (rg->idx < i) {
    esmac_rng_next(rg);
  }
  return rg->prn;
}

// uniform in [a,b)
double esmac_rng_get_double(esmac_rng_t *rg, double a, double b, int64_t i) {
  uint64_t ui = esmac_rng_get(rg,i);
  unsigned long int u1 = ui & 0xFFFFFFFF;
  unsigned long int u2 = (ui >> 32) & 0xFFFFFFFF;
  double v = ldexp((double) u1,-64) + ldexp((double) u2,-32);
  return a + (b-a)*v;
}

uint64_t esmac_generate_seed(int short_seed) {
  int ti = (int)time(NULL);  // do whatever you want here
  return (uint64_t) short_seed + ( ((uint64_t) ti) << 32 );
}
