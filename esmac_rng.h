#ifndef ESMAC_RNG_H
#define ESMAC_RNG_H

#include <stdint.h>

typedef struct {
  uint64_t s[16];
  int p;
  
  uint64_t seed;
  uint64_t sinit[16];
  
  uint64_t prn;
  int64_t idx;
  
} esmac_rng_t;

esmac_rng_t * esmac_rng_alloc(uint64_t seed);
void esmac_rng_free(esmac_rng_t * rg);

void esmac_rng_reset(uint64_t seed, esmac_rng_t *rg);

// xorshift random number generator according to Sebastiano Vigna, arxiv:1402.6246
// [ ACM Trans. Math. Software 42, 4 (2016) ]
uint64_t esmac_rng_next(esmac_rng_t *rg);

uint64_t esmac_rng_get(esmac_rng_t * rg, int64_t i);

// uniform in [a,b)
double esmac_rng_get_double(esmac_rng_t *rg, double a, double b, int64_t i);

/*
// rv[i] stores for i0 <= i < i1
// length rv >= i1-i0
int esmac_ranvec_int(uint64_t seed, int i0, int i1, uint64_t *rv);
int esmac_ranvec_double(uint64_t seed, int i0, int i1, double *rv);
*/

uint64_t esmac_generate_seed(int short_seed);

#endif /* ESMAC_RNG_H */
