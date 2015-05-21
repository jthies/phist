#include "phist_config.h"
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "gtest/gtest.h"

#include "phist_macros.h"
#include <cmath>


#if defined(PHIST_KERNEL_LIB_BUILTIN) && defined(PHIST_HIGH_PRECISION_KERNELS)
#include "prec_helpers.h"


TEST(BuiltinPrecTest, DOUBLE_2SUM)
{
  double a, b, s, t;

  a = 1, b = 0;
  DOUBLE_2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(0,t);

  a = 99; b = 523;
  DOUBLE_2SUM(a,b,s,t);
  EXPECT_EQ(a+b,s);
  EXPECT_EQ(0,t);

  a = 1, b = 1.e-32;
  DOUBLE_2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(b,t);

  a = 1.e-32, b = 1;
  DOUBLE_2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(a,t);

  a = (1<<8), b = -(1<<8) + 1;
  DOUBLE_2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(0,t);
}


TEST(BuiltinPrecTest, DOUBLE_FAST2SUM)
{
  double a, b, s, t;

  a = 1, b = 0;
  DOUBLE_FAST2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(0,t);

  a = 99; b = 523;
  DOUBLE_2SUM(a,b,s,t);
  EXPECT_EQ(a+b,s);
  EXPECT_EQ(0,t);

  a = 1, b = 1.e-32;
  DOUBLE_FAST2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(b,t);

  a = (1<<8), b = -(1<<8) + 1;
  DOUBLE_FAST2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(0,t);

  a = 1.e-32, b = 1;
  DOUBLE_FAST2SUM(a,b,s,t);
  EXPECT_EQ(1,s);
  EXPECT_EQ(0,t); // should be 'a', but FAST2SUM is not accurate enough
}


TEST(BuiltinPrecTest, DOUBLE_2MULTFMA)
{
  double a, b, s, t;

  a = 37; b = -92;
  DOUBLE_2MULTFMA(a,b,s,t);
  EXPECT_EQ(a*b,s);
  EXPECT_EQ(0,t);

  a = 1./(1l<<50l)+1, b = 1./(1l<<50l)+1;
  DOUBLE_2MULTFMA(a,b,s,t);
  //PHIST_SOUT(PHIST_INFO,"%e * %e = %e + %e\n", a, b, s, t);
  EXPECT_EQ(2./(1l<<50l)+1,s);
  EXPECT_EQ(1./(1l<<50l)*1./(1l<<50l),t);
}


TEST(BuiltinPrecTest, DOUBLE_4SUM)
{
  double a, aC, b, bC, s, t;

  a = 1, aC = 0, b = 1, bC = 0;
  DOUBLE_4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(2,s);
  EXPECT_EQ(0,t);

  a = (1<<8), aC = -1, b = (1<<10), bC = 2;
  DOUBLE_4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(a+aC+b+bC,s);
  EXPECT_EQ(0,t);
  DOUBLE_4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(a+aC+b+bC,s);
  EXPECT_EQ(0,t);

  a = (1l<<60), aC = 1., b = 1., bC = 1./(1l<<60);
  DOUBLE_4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(a,s);
  EXPECT_EQ(aC+b,t);
  DOUBLE_4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(a,s);
  EXPECT_EQ(aC+b,t);

  a = -(1l<<60), aC = 2., b = (1l<<60), bC = 1.;
  DOUBLE_4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(aC+bC,s);
  EXPECT_EQ(0,t);
  DOUBLE_4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(aC+bC,s);
  EXPECT_EQ(0,t);

  a = -(1l<<60), aC = 2., b = (1l<<60), bC = 1./(1l<<60);
  DOUBLE_4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(aC,s);
  EXPECT_EQ(bC,t);
  DOUBLE_4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(aC,s);
  EXPECT_EQ(bC,t);
}


TEST(BuiltinPrecTest, DOUBLE_FAST4SUM)
{
  double a, aC, b, bC, s, t;

  a = 1, aC = 0, b = 1, bC = 0;
  DOUBLE_FAST4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(2,s);
  EXPECT_EQ(0,t);

  a = (1<<8), aC = -1, b = (1<<10), bC = 2;
  DOUBLE_FAST4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(a+aC+b+bC,s);
  EXPECT_EQ(0,t);
  DOUBLE_FAST4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(a+aC+b+bC,s);
  EXPECT_EQ(0,t);

  a = (1l<<60), aC = 1., b = 1., bC = 1./(1l<<60);
  DOUBLE_FAST4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(a,s);
  EXPECT_EQ(aC+b,t);
  DOUBLE_FAST4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(a,s);
  EXPECT_EQ(aC,t); // should be aC+b, but not accurate enough

  a = -(1l<<60), aC = 2., b = (1l<<60), bC = 1.;
  DOUBLE_FAST4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(aC+bC,s);
  EXPECT_EQ(0,t);
  DOUBLE_FAST4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(aC+bC,s);
  EXPECT_EQ(0,t);

  a = -(1l<<60), aC = 2., b = (1l<<60), bC = 1./(1l<<60);
  DOUBLE_FAST4SUM(a,aC,b,bC,s,t);
  EXPECT_EQ(aC,s);
  EXPECT_EQ(bC,t);
  DOUBLE_FAST4SUM(b,bC,a,aC,s,t);
  EXPECT_EQ(aC,s);
  EXPECT_EQ(0,t); // should be bC, but not accurate enough
}


TEST(BuiltinPrecTest, DOUBLE_4MULTFMA)
{
  double a, aC, b, bC, s, t;

  a = 99, aC = 1, b = 7, bC = 1;
  DOUBLE_4MULTFMA(a,aC,b,bC,s,t);
  EXPECT_EQ(800,s+t);

  a = (1l<<30), aC = 1./(1l<<30), b = (1l<<30), bC = 1./(1l<<30);
  DOUBLE_4MULTFMA(a,aC,b,bC,s,t);
  EXPECT_EQ((1l<<60),s);
  EXPECT_EQ(2,t);
}


TEST(BuiltinPrecTest, DOUBLE_2DIVFMA)
{
  double a, b, s, t, s_, t_;

  a = -2468, b = 2;
  DOUBLE_2DIVFMA(a,b,s,t);
  EXPECT_EQ(-1234,s);
  EXPECT_EQ(0,t);

  a = 1., b = 3.;
  DOUBLE_2DIVFMA(a,b,s,t);
  PHIST_SOUT(PHIST_INFO,"%e / %e = %e + %e\n", a, b, s, t);
  EXPECT_EQ(1./3.,s);
  DOUBLE_4MULTFMA(s,t,b,0.,s_,t_);
  PHIST_SOUT(PHIST_INFO,"( %e + %e ) * %e = %e + %e\n", s, t, b, s_, t_);
  EXPECT_EQ(1.,s_);
  EXPECT_LE(std::abs(t_), 1.e-30);

}


TEST(BuiltinPrecTest, DOUBLE_4DIV_NEWTONRAPHSON_FMA)
{
  double a, aC, s, t, s_, t_;

  a = 128, aC = 0;
  DOUBLE_4DIV_NEWTONRAPHSON_FMA(a,aC,s,t);
  PHIST_SOUT(PHIST_INFO," 1 / ( %e + %e ) = %e + %e\n", a, aC, s, t);
  EXPECT_EQ(1./128,s);
  EXPECT_LE(std::abs(t), 1.e-30);

  a = 127, aC = 1;
  DOUBLE_4DIV_NEWTONRAPHSON_FMA(a,aC,s,t);
  PHIST_SOUT(PHIST_INFO," 1 / ( %e + %e ) = %e + %e\n", a, aC, s, t);
  EXPECT_EQ(1./128,s);
  EXPECT_LE(std::abs(t), 1.e-20);

  a = 3, aC = 0;
  DOUBLE_4DIV_NEWTONRAPHSON_FMA(a,aC,s,t);
  EXPECT_EQ(1./3.,s);
  DOUBLE_4MULTFMA(s,t,a,aC,s_,t_);
  PHIST_SOUT(PHIST_INFO,"( %e + %e ) * ( %e + %e ) = %e + %e\n", s, t, a, aC, s_, t_);
  EXPECT_EQ(1.,s_);
  EXPECT_LE(std::abs(t_), 1.e-20);

  DOUBLE_2DIVFMA(1.,3.,a,aC);
  DOUBLE_4DIV_NEWTONRAPHSON_FMA(a,aC,s,t);
  PHIST_SOUT(PHIST_INFO," 1 / ( %e + %e ) = %e + %e\n", a, aC, s, t);
  EXPECT_EQ(3.,s);
  EXPECT_LE(std::abs(t), 1.e-30);
}


TEST(BuiltinPrecTest, DOUBLE_2SQRTFMA)
{
  double a, s, t, s_, t_;

  a = 4;
  DOUBLE_2SQRTFMA(a,s,t);
  EXPECT_EQ(2,s);
  EXPECT_EQ(0,t);

  a = 2;
  DOUBLE_2SQRTFMA(a,s,t);
  PHIST_SOUT(PHIST_INFO,"sqrt(%e) = %e + %e\n", a, s, t);
  EXPECT_EQ(sqrt(2.),s);
  DOUBLE_4MULTFMA(s,t,s,t,s_,t_);
  PHIST_SOUT(PHIST_INFO,"(%e + %e) * (%e + %e) = %e + %e\n", s, t, s, t, s_, t_);
  EXPECT_EQ(2.,s_);
  EXPECT_LE(std::abs(t_),1.e-30);
}


TEST(BuiltinPrecTest, DOUBLE_4SQRT_NEWTONRAPHSON_FMA)
{
  double a, aC, sqrt_a, sqrt_aC, divsqrt_a, divsqrt_aC, s, t;

  a = 16, aC = 0;
  DOUBLE_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC);
  EXPECT_EQ(4,sqrt_a);
  EXPECT_EQ(0,sqrt_aC);
  EXPECT_EQ(0.25,divsqrt_a);
  EXPECT_EQ(0,divsqrt_aC);

  a = 255, aC = 1;
  DOUBLE_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC);
  PHIST_SOUT(PHIST_INFO,"sqrt(%e+%e) = %e + %e\n", a, aC, sqrt_a, sqrt_aC);
  PHIST_SOUT(PHIST_INFO,"1/sqrt(%e+%e) = %e + %e\n", a, aC, divsqrt_a, divsqrt_aC);
  EXPECT_EQ(16,sqrt_a);
  EXPECT_LE(std::abs(sqrt_aC),1.e-16);
  EXPECT_EQ(1./16,divsqrt_a);
  EXPECT_LE(std::abs(divsqrt_aC),1.e-16);

  a = (1l<<30)+1, aC = 1./(1l<<30)+1./(1l<<60);
  DOUBLE_4SQRT_NEWTONRAPHSON_FMA(a,aC,sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC);
  PHIST_SOUT(PHIST_INFO,"sqrt(%e+%e) = %e + %e\n", a, aC, sqrt_a, sqrt_aC);
  PHIST_SOUT(PHIST_INFO,"1/sqrt(%e+%e) = %e + %e\n", a, aC, divsqrt_a, divsqrt_aC);
  DOUBLE_4MULTFMA(sqrt_a,sqrt_aC,sqrt_a,sqrt_aC,s,t);
  PHIST_SOUT(PHIST_INFO,"(%e + %e) * (%e + %e) = %e + %e\n", sqrt_a, sqrt_aC, sqrt_a, sqrt_aC, s, t);
  EXPECT_EQ(a,s);
  EXPECT_EQ(aC,t);
  DOUBLE_4MULTFMA(sqrt_a,sqrt_aC,divsqrt_a,divsqrt_aC,s,t);
  PHIST_SOUT(PHIST_INFO,"(%e + %e) * (%e + %e) = %e + %e\n", sqrt_a, sqrt_aC, divsqrt_a, divsqrt_aC, s, t);
  EXPECT_EQ(1,s);
  EXPECT_LE(std::abs(t),1.e-30);
}


#endif /* PHIST_KERNEL_LIB BUILTIN && PHIST_HIGH_PRECISION_KERNELS */


