#ifndef PHIST_CORE_FLAGS_H
#define PHIST_CORE_FLAGS_H

/* some flags that influence the behavior of core functions, see also
   phist_kernel_flags.h for an introduction to flags in PHIST.
   
   The values 2^[16...23] are reserved for core functionality for now, again, 
   if there is no function that accepts flags A and B, A==B is allowed.
   
 */

/* for KPM */
#define PHIST_KPM_SINGLEVEC 65536

/* for orthog routines: fill up output vector with random numbers and orthogonalize them along */
/* with the original entries if the vector W is found to be rank-deficient.                    */
#define PHIST_ORTHOG_RANDOMIZE_NULLSPACE 131072

#endif
