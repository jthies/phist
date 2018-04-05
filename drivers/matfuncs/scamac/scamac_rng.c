#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>

#include "scamac_rng.h"

#include "scamac_rng_jump1024table_inc.c"

/* ++++++++++ internal functions +++++++++++ */

// generate 1024 pseudo-random values, starting from state st, and store into ranv.
// st is replaced by "next^1024 st" after return.
static void generate_1024(uint64_t * st, uint64_t *ranv) {
  int i;
  for (i=0;i<1024;i++) {
    ranv[i]=st[i%16] * UINT64_C(1181783497276652981); // scramble
    // next
    uint64_t s0 = st[i%16];
    uint64_t s1 = st[(i+1)%16];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    st[(i+1)%16] =  s0 ^ s1;
  }
}

// jump ahead by 2^n steps, using precomputed table values in jump1024.txt
// Powers n<10 are trivial, since smaller than the dimension 1024 of the xorshift generator,
// and could also be accomplished by calling next() the respective number of times. 
static void jump_ahead_pow2(int n,  uint64_t * st) {
  assert( (n>=0) && (n<=512) );
  int i,j,k;
  uint64_t tt[16];
  for (k=0;k<16;k++) {
    tt[k]=st[k];
    st[k]=0;
  }
  for (i=0;i<16;i++) {
    for (j=0;j<64;j++) {
	    if (JUMP1024[16*n+i] & (UINT64_C(1) << j)) {
	      for (k=0;k<16;k++) {
         st[k] ^= tt[(k+j)%16]; // check
        //  st[(k+j)%16] ^= tt[k]; // check
	      }
	    }
	  // next
	  uint64_t s0 = tt[j%16];
	  uint64_t s1 = tt[(j+1)%16];
	  s1 ^= s1 << 31; // a
	  s1 ^= s1 >> 11; // b
	  s0 ^= s0 >> 30; // c
	  tt[(j+1)%16] =  s0 ^ s1;
    }
  }
}

/*
// jump ahead by one step.
// Mainly used for testing, includes some unnecessary copying. 
static void next(uint64_t * st) {
  uint64_t s0 = st[0];
  uint64_t s1 = st[1];
  s1 ^= s1 << 31; // a
  s1 ^= s1 >> 11; // b 
  s0 ^= s0 >> 30; // c
  st[1] = s0 ^ s1;
  // and swap
  s0 = st[0];
  st[ 0]=st[ 1];
  st[ 1]=st[ 2];
  st[ 2]=st[ 3];
  st[ 3]=st[ 4];
  st[ 4]=st[ 5];
  st[ 5]=st[ 6];
  st[ 6]=st[ 7];
  st[ 7]=st[ 8];
  st[ 8]=st[ 9];
  st[ 9]=st[10];
  st[10]=st[11];
  st[11]=st[12];
  st[12]=st[13];
  st[13]=st[14];
  st[14]=st[15];
  st[15]=s0;
}
*/

static void jump_ahead(uint64_t n,  uint64_t * st) {
  int i=0;
  while (n) {
    if (n & 1) {
      jump_ahead_pow2(i, st);
    }
    i++;
    n = n >> 1;
  }
}

// seed with a uint64 seed
static void seed_to_state(uint64_t seed, uint64_t * st) {
  int i;
  for (i=0;i<16;i++) {
    st[i]=seed;
    seed=(~seed) * UINT64_C(12530480023924335475); // just some number
  }
  // jump ahead by 2^37+2^27+2^17 (just some number)
  jump_ahead_pow2(17, st);
  jump_ahead_pow2(27, st);
  jump_ahead_pow2(37, st);
}

/* +++++++++++++++++++++++++++++++++ */

ScamacErrorCode scamac_ransrc_alloc(scamac_rng_seed_ty seed, scamac_ransrc_st ** ransrc) {
  if (!ransrc) { return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT; }
  scamac_ransrc_st * rs = malloc(sizeof *rs);
  if (!rs) { return SCAMAC_EMALLOCFAIL; }
  ScamacErrorCode err;
  err = scamac_ransrc_set_seed(seed, rs);
  if (err) {
    free(rs);
    return err;
  }
  *ransrc = rs;
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_ransrc_free(scamac_ransrc_st * ransrc) {
//  if (ransrc) {
    free(ransrc);
//  }
    return SCAMAC_EOK;
}

ScamacErrorCode scamac_ransrc_set_seed(scamac_rng_seed_ty seed, scamac_ransrc_st * ransrc) {
  if (!ransrc) { return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT; }
  seed_to_state(seed, ransrc->sinit);
  ransrc->sliceidx=-1;
  return SCAMAC_EOK;
}

// xorshift random number generator according to Sebastiano Vigna, arxiv:1402.6246
// [ ACM Trans. Math. Software 42, 4 (2016) ]
/*
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
*/

uint64_t scamac_ransrc_uint64(scamac_ransrc_st * ransrc, uint64_t i) {
  // in bunches of 1024 numbers
  if ( (i >> 10) != ransrc->sliceidx ) {//generate first
    uint64_t st[16];
    int j;
    for (j=0;j<16;j++) { st[j]=ransrc->sinit[j]; }
    //  jump_ahead(i >> 10, st);
    jump_ahead( (i >> 10) << 10, st); // to multiple of 1024
    generate_1024(st, ransrc->prns);
    ransrc->sliceidx = i >> 10;
  }
  return ransrc->prns[i % 1024];
}

// uniform in [a,b)
double scamac_ransrc_double(scamac_ransrc_st *ransrc, double a, double b, uint64_t i) {
  uint64_t u = scamac_ransrc_uint64(ransrc,i);
  // uint64 has 64 bits of information, but double has only 53 mantissa bits.
  // For a uniformly distributed number, we lose the 11 bits.
  double v = ( (double) (u>>11) ) / (1.0 + (UINT64_MAX >> 11));
  return a + (b-a)*v;
}


ScamacErrorCode scamac_ranvec_uint64(scamac_rng_seed_ty seed, uint64_t iv, uint64_t nv, uint64_t * ranvec) {
  if (!ranvec) { return SCAMAC_ENULL | 4 << SCAMAC_ESHIFT; }

  uint64_t st[16];
  seed_to_state(seed, st);

  // multiple of 1024 below i

  uint64_t mi;
  mi = (iv >> 10) << 10;

  jump_ahead(mi, st);

  uint64_t i,j;

  // skip first
  for (i=mi;i<iv;i++) {
    // next
    uint64_t s0 = st[i%16];
    uint64_t s1 = st[(i+1)%16];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    st[(i+1)%16] =  s0 ^ s1;
  }

  // now, produce
  j=0;
  for (i=iv;i<iv+nv;i++) {
    ranvec[j]=st[i%16] * UINT64_C(1181783497276652981); // scramble
    j++;
    // next
    uint64_t s0 = st[i%16];
    uint64_t s1 = st[(i+1)%16];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    st[(i+1)%16] =  s0 ^ s1;
  }  

  return SCAMAC_EOK;

}

ScamacErrorCode scamac_ranvec_double(scamac_rng_seed_ty seed, uint64_t iv, uint64_t nv, double a, double b, double * ranvec) {
  if (!ranvec) { return SCAMAC_ENULL | 4 << SCAMAC_ESHIFT; }

  uint64_t st[16];
  seed_to_state(seed, st);

  // multiple of 1024 below i

  uint64_t mi;
  mi = (iv >> 10) << 10;

  jump_ahead(mi, st);

  uint64_t i,j;

  // skip first
  for (i=mi;i<iv;i++) {
    // next
    uint64_t s0 = st[i%16];
    uint64_t s1 = st[(i+1)%16];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    st[(i+1)%16] =  s0 ^ s1;
  }

  // now, produce
  j=0;
  for (i=iv;i<iv+nv;i++) {
    uint64_t u = st[i%16] * UINT64_C(1181783497276652981); // scramble
    double v = ( (double) (u>>11) ) / (1.0 + (UINT64_MAX >> 11));
    ranvec[j] = a + (b-a)*v;
    j++;
    // next
    uint64_t s0 = st[i%16];
    uint64_t s1 = st[(i+1)%16];
    s1 ^= s1 << 31; // a
    s1 ^= s1 >> 11; // b
    s0 ^= s0 >> 30; // c
    st[(i+1)%16] =  s0 ^ s1;
  }  

  return SCAMAC_EOK;

}



scamac_rng_seed_ty scamac_rng_string_to_seed(const char *str) {
  size_t l = strlen(str);
  // skip initial blanks
  while (l>0) {
    if (str[0] == ' ') {
      str++;
      l-- ;
    } else {
      break;
    }
  }
  if (l<1) { // fallback, really
    return 0;
  }
  if (l<=20) { // could be a decimal number of uint64
    int is_dec=1;
    int i;
    for (i=0;i<l;i++) {
      if (str[i]<'0' || str[i]>'9') {
        is_dec=0;
        break;
      }
    }
    if (is_dec) {
      uint64_t seed;
      sscanf(str,"%"SCNu64,&seed);
      return seed;
    }
  }
  if (l<=18) { // could be a hex number of uint64
    int is_hex=1;
    if (l>2) {// skip initial "0x", if present
      if (str[0]=='0' && (str[1]=='x' || str[1]=='X')) {
        str += 2;
        l -= 2;
      }
    }
    int i;
    for (i=0;i<l;i++) {
      if (str[i]<'0' || str[i]>'f') {
        is_hex=0;
        break;
      }
      if (str[i]>'9' && str[i]<'A') {
        is_hex=0;
        break;
      }
      if (str[i]>'F' && str[i]<'a') {
        is_hex=0;
        break;
      }
    }
    if (is_hex) {
      uint64_t seed;
      sscanf(str,"%"SCNx64,&seed);
      return seed;
    }
  }
  
  uint64_t seed=0;
  int i;
  for (i=0;i<l;i++) {
    seed= ((seed << 11) | (seed >> 53)) ^ str[i];
  }
  return seed;
}
