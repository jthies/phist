/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifdef PHISTTEST_MVEC_CREATE
#undef PHISTTEST_MVEC_CREATE
#endif

#ifdef PHISTTEST_MVEC_CLONE_SHAPE
#undef PHISTTEST_MVEC_CLONE_SHAPE
#endif

#ifdef PHISTTEST_RMVEC_CREATE
#undef PHISTTEST_RMVEC_CREATE
#endif

#ifdef PHISTTEST_CMVEC_CREATE
#undef PHISTTEST_CMVEC_CREATE
#endif

//! macro to create an mvec in the currently defined numerical type,
//! pass in correct flags for testing (e.g. replicate device memory on
//! CPU)
#define PHISTTEST_MVEC_CREATE(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  SUBR(mvec_create)(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR);

#define PHISTTEST_MVEC_CLONE_SHAPE(MVEC_STARSTAR,MVEC_STAR,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  SUBR(mvec_clone_shape)(MVEC_STARSTAR,MVEC_STAR,IFLAG_STAR);

#ifdef IS_DOUBLE

# define PHISTTEST_RMVEC_CREATE(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  phist_Dmvec_create(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR);

# define PHISTTEST_CMVEC_CREATE(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  phist_Zmvec_create(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR);

#else

# define PHISTTEST_RMVEC_CREATE(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  phist_Smvec_create(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR);

# define PHISTTEST_CMVEC_CREATE(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  phist_Cmvec_create(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR);

#endif
