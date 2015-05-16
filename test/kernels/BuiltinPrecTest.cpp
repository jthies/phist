#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"


#if defined(PHIST_KERNEL_LIB_BUILTIN) && defined(PHIST_HIGH_PRECISION_KERNELS)
#include "prec_helpers.h"


TEST(BuiltinPrecTest, DOUBLE_2SUM)
{
  double a, b, s, t;

  a = 1, b = 0;
  DOUBLE_2SUM(a,b,s,t);
  ASSERT_EQ(s,a);
  ASSERT_EQ(t,b);


  a = 1, b = 1.e-32;
  DOUBLE_2SUM(a,b,s,t);
  ASSERT_EQ(s,a);
  ASSERT_EQ(t,b);

  a = 1.e-32, b = 1;
  DOUBLE_2SUM(a,b,s,t);
  ASSERT_EQ(s,b);
  ASSERT_EQ(t,a);
}

#endif /* PHIST_KERNEL_LIB BUILTIN && PHIST_HIGH_PRECISION_KERNELS */


