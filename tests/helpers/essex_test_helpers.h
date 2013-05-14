#ifndef ESSEX_TEST_HELPERS_H
#define ESSEX_TEST_HELPERS_H

//! for testing type safety
typedef struct
  {
  int i_component;
  double d_component;
  float* fp_component;
  } essex_bad_object;

#endif
