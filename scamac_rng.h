/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_RNG_H
#define SCAMAC_RNG_H

#include <stdint.h>

typedef struct {
  uint64_t s[16];
  int p;

  uint64_t seed;
  uint64_t sinit[16];

  uint64_t prn;
  int64_t idx;

} scamac_rng_st;

scamac_rng_st * scamac_rng_alloc(uint64_t seed);
void scamac_rng_free(scamac_rng_st * rg);

void scamac_rng_reset(uint64_t seed, scamac_rng_st *rg);

// xorshift random number generator according to Sebastiano Vigna, arxiv:1402.6246
// [ ACM Trans. Math. Software 42, 4 (2016) ]
uint64_t scamac_rng_next(scamac_rng_st *rg);

uint64_t scamac_rng_get(scamac_rng_st * rg, int64_t i);

// uniform in [a,b)
double scamac_rng_get_double(scamac_rng_st *rg, double a, double b, int64_t i);

/*
// rv[i] stores for i0 <= i < i1
// length rv >= i1-i0
int scamac_ranvec_int(uint64_t seed, int i0, int i1, uint64_t *rv);
int scamac_ranvec_double(uint64_t seed, int i0, int i1, double *rv);
*/

uint64_t scamac_generate_seed(int short_seed);

#endif /* SCAMAC_RNG_H */
