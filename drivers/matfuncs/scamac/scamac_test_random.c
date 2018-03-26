#include <stdio.h>
#include <inttypes.h>

#include <scamac_rng.h>


int main(int argc, char * argv[]) {

  scamac_rng_seed_ty seed;

  // some arbitary string as seed
  seed = scamac_rng_string_to_seed("uihehr");

  scamac_ransrc_st * ransrc;
  scamac_ransrc_alloc(seed, &ransrc);

  uint64_t i1,i2,n;

  n=1000;
  i1=1234;
  i2=i1+n;

  uint64_t * rv1, *rv2;
  rv1 = malloc(2 * n * sizeof * rv1);
  rv2 = malloc(2 * n * sizeof * rv2);

  scamac_ranvec_uint64(seed, i1, 2*n, rv1);
  scamac_ranvec_uint64(seed, i2, 2*n, rv2);

 uint64_t i;

 for (i=0;i<n;i++) {    
    if (rv1[i+n] != rv2[i] || rv1[i+n] != scamac_ransrc_uint64(ransrc, i1+i+n)) {
      printf("Not identical. ERROR.\n");
      return 1;
    }
  }

  printf("Identical. Success.\n");

  scamac_ransrc_free(ransrc);
  ransrc=NULL;

  return 0;
}
