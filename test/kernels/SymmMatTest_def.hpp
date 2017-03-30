/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME>,
                 public KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,2> 
{

public:
  typedef KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,2> VTest;

  static void SetUpTestCase()
  {
    int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
    SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
    VTest::SetUpTestCase();
  }

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    SparseMatTest::SetUp();
    VTest::SetUp();

    haveMat_ = false;
    if (typeImplemented_ && !problemTooSmall_)
    {
      haveMat_ = (A_ != NULL);
    }
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    VTest::TearDown();
    SparseMatTest::TearDown();
  }

  static void TearDownTestCase()
  {
    VTest::TearDownTestCase();
    SparseMatTest::TearDownTestCase();
  }

protected:

  bool haveMat_;
};

  TEST_F(CLASSNAME, read_matrices) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_TRUE(AssertNotNull(A_));
    }
  }

  TEST_F(CLASSNAME, A1_probe_symmetry)
  {
    if (typeImplemented_ && !problemTooSmall_ && haveMat_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
      TYPE(sdMat_ptr) VtAV=NULL;
      SUBR(sdMat_create)(&VtAV,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),vec1_,vec2_,st::zero(),VtAV,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(VtAV,&iflag_);
#endif
      ASSERT_EQ(0,iflag_);

      _ST_* C_raw=NULL;
      int lda;
      SUBR(sdMat_extract_view)(VtAV,&C_raw,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check V'AV != 0 and V'AV symmetric      
      _MT_ sum=mt::one();;
      _MT_ diff=mt::zero();
      for (int i=0;i<_NV_;i++)
        for (int j=i;j<_NV_;j++)
        {
          _ST_ cij=C_raw[i*lda+j];
          _ST_ cji=C_raw[j*lda+i];
          sum=std::min(sum,st::abs(cij)+st::abs(cji));
          diff=std::max(diff,st::abs(st::conj(cij)-cji));
        }
      ASSERT_NEAR(mt::zero(),diff,sqrt(mt::eps()));
      ASSERT_TRUE(sum>sqrt(mt::eps()));

      SUBR(sdMat_delete)(VtAV,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
