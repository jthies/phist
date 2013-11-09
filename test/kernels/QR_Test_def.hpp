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
      SUBR(mvec_random)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<stride_*nloc_;i+=stride_)
          {
          vec2_vp_[j*lda_+i] = vec1_vp_[j*lda_+i];
          }
        }
//      PrintVector(*cout,"QR_Test V",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      SUBR(mvec_QR)(vec2_,mat1_,&ierr_);
//      PrintVector(*cout,"QR_Test Q",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
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
    if (typeImplemented_)
      {
      ASSERT_REAL_EQ(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_));
      }
    }


  // check if vectors are mutually orthogonal after QR factorization
  TEST_F(CLASSNAME, ortho) 
    {
    if (typeImplemented_)
      {
      ASSERT_REAL_EQ(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_));
      }
    }

  TEST_F(CLASSNAME, with_rank_deficiency) 
    {
    if (typeImplemented_)
      {
      if (nvec_==1)
        {
        SUBR(mvec_put_value)(vec1_,st::zero(),&ierr_);
      ASSERT_EQ(0,ierr_);
        }
      else
        {
        SUBR(mvec_random)(vec1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        // set last two columns to same vector
        for (int i=0;i<stride_*nloc_;i+=stride_)
          {
          vec1_vp_[(nvec_-1)*lda_+i] = vec1_vp_[(nvec_-2)*lda_+i];
          }
        }
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_QR)(vec2_,mat1_,&ierr_);
      // check that the rank deficiency was detected
      ASSERT_EQ(1, ierr_);
      // check that we anyway got something orthogonal back
      ASSERT_REAL_EQ(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_));
      ASSERT_REAL_EQ(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_));
      }
    }

