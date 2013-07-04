#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_> 
  {

public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::SetUp();
    if (typeImplemented_)
      {
      _SUBR_(mvec_random)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<stride_*nrows_;i+=stride_)
          {
          vec2_vp_[j*lda_+i] = vec1_vp_[j*lda_+i];
          }
        }
      _SUBR_(mvec_QR)(vec2_,mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::TearDown();
    }

};

  // check if vectors are normalized correctly after QR factorization
  TEST_F(CLASSNAME, normality) 
    {
    ASSERT_REAL_EQ((_MT_)1.0,ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_));
    }


  // check if vectors are mutually orthogonal after QR factorization
  TEST_F(CLASSNAME, ortho) 
    {
    ASSERT_REAL_EQ((_MT_)1.0,ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_));
    }

