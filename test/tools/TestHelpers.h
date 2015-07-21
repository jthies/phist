#ifdef PHISTTEST_MVEC_CREATE
#undef PHISTTEST_MVEC_CREATE
#endif

//! macro to create an mvec in the currently defined numerical type,
//! pass in correct flags for testing (e.g. replicate device memory on
//! CPU)
#define PHISTTEST_MVEC_CREATE(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR) \
  *(IFLAG_STAR)=PHIST_MVEC_REPLICATE_DEVICE_MEM; \
  SUBR(mvec_create)(MVEC_STARSTAR,MAP_STAR,NVEC_IN,IFLAG_STAR); \

