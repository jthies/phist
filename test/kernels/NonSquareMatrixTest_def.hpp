/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

using namespace phist::testing;

/*
#ifndef CONCAT
#define CONCAT(a,b) a ## b
#endif
#ifdef MTEST_NAME
#undef MTEST_NAME
#endif
#define VMTEST CONCAT(CLASSNAME,_VmTest) 

class VMTEST
: public KernelTestWithVectors<_ST_,_M_,_NV_,_USE_VIEWS_,3>
{
  VMTEST()
  {
  }
  ~VMTEST()
  {
  }
};
*/

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,_M_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,3>,
                 public virtual KernelTestWithVectors<_ST_,_M_,_NV_,_USE_VIEWS_,3>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_,_USE_VIEWS_>
{

  public:
  
  typedef KernelTestWithSparseMat<_ST_,_N_,_M_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,3> VnTest;
  typedef KernelTestWithVectors<_ST_,_M_,_NV_,_USE_VIEWS_,3> VmTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_,_USE_VIEWS_> MTest;
  typedef TestWithType< _MT_ > MT_Test;

  static void SetUpTestCase()
  {
    int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
    SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
    VnTest::SetUpTestCase();
    VmTest::SetUpTestCase();
    MTest::SetUpTestCase();
  }

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    SparseMatTest::SetUp();
    VnTest::SetUp();
    VmTest::SetUp();
    MTest::SetUp();
    
    haveMat_ = (A_ != NULL);
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    MTest::TearDown();
    VnTest::TearDown();
    VmTest::TearDown();
    SparseMatTest::TearDown();
  }

  static void TearDownTestCase()
  {
    MTest::TearDownTestCase();
    VnTest::TearDownTestCase();
    VmTest::TearDownTestCase();
    SparseMatTest::TearDownTestCase();
  }

protected:

  bool haveMat_;
};

  TEST_F(CLASSNAME, read_matrices) 
  {
    if (typeImplemented_ && !VnTest::problemTooSmall_)
    {
      ASSERT_TRUE(AssertNotNull(A_));
    
      // test that the global number of rows/cols is correct in the objects
      phist_gidx gnrows, gncols;
      SUBR(sparseMat_global_nrows)(A_,&gnrows,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(gnrows,VnTest::nglob_);
      SUBR(sparseMat_global_ncols)(A_,&gncols,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(gncols,VmTest::nglob_);
    }
  }


#if MATNAME == MATNAME_IDFUNC
  TEST_F(CLASSNAME, I_pad_times_mvec)
  {
    // the matrix is the identity matrix padded by 0, e.g. [eye(n) zeros(n,m-n)] or [I;zeros(n-m,m)]
    // if m>n or n>m, respectively.
    if (!typeImplemented_ || VnTest::problemTooSmall_ || !haveMat_) return;
    
    ST alpha, beta;
    //I*X=X?
    SUBR(mvec_random)(VnTest::vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(VmTest::vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),VnTest::vec1_,st::zero(),VnTest::vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),VmTest::vec1_,st::zero(),VmTest::vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(sparseMat_times_mvec)(st::one(),A_,VmTest::vec1_,st::zero(),VnTest::vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    // Testing the results is non-trivial because the parallel distribution must be taken into account.
    // The first min(n,m) rows of the two vectors should the same.
    // ...
#if _NROWS_>_NCOLS_
    // The last n-m rows of VnTest::vec1_ must be 0.
    // ...
#endif
  }
#endif
