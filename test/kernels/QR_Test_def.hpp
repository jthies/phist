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
      ASSERTREALEQ(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_));
      }
    }


  // check if vectors are mutually orthogonal after QR factorization
  TEST_F(CLASSNAME, ortho) 
    {
    if (typeImplemented_)
      {
      ASSERTREALEQ(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_));
      }
    }

