#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_> 
  {

public:

  typedef KernelTestWithVectors<_ST_,_N_,_NV_> VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_> MTest;

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    VTest::SetUp();
    MTest::SetUp();
    if (typeImplemented_)
      {
      SUBR(mvec_random)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<stride_*nloc_;i+=stride_)
          {
#ifdef PHIST_MVECS_ROW_MAJOR
          vec2_vp_[j+i*lda_] = vec1_vp_[j+i*lda_];
#else
          vec2_vp_[j*lda_+i] = vec1_vp_[j*lda_+i];
#endif
          }
        }
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    VTest::TearDown();
    MTest::TearDown();
    }

  _MT_ nrms_[_NV_];
};

  TEST_F(CLASSNAME, DISABLED_with_random_vectors) 
    {
    if (typeImplemented_)
      {
//      PrintVector(*cout,"SVQB_Test V",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      SUBR(svqb)(vec2_,mat1_,nrms_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // doing a SVQB decomp must not relocate data:
      ASSERT_EQ(true,MTest::pointerUnchanged(mat1_,mat1_vp_,m_lda_));
      ASSERT_EQ(true,VTest::pointerUnchanged(vec2_,vec2_vp_,lda_));
//      PrintVector(*cout,"SVQB_Test Q",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      // check norms
      _MT_ nrms_ref[_NV_];
      SUBR(mvec_norm2)(vec1_,nrms_ref,&ierr_);
      ASSERT_EQ(0,ierr_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_REAL_EQ(nrms_ref[i], nrms_[i]);
      }
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));

      // check Q=V*B
      SUBR(mvec_times_sdMat)(-st::one(),vec1_,mat1_,st::one(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
      }
    }

  TEST_F(CLASSNAME, DISABLED_with_rank_deficiency) 
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
#ifdef PHIST_MVECS_ROW_MAJOR
          vec1_vp_[(nvec_-1)+i*lda_] = vec1_vp_[(nvec_-2)+i*lda_];
#else
          vec1_vp_[(nvec_-1)*lda_+i] = vec1_vp_[(nvec_-2)*lda_+i];
#endif
          }
        }
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(svqb)(vec2_,mat1_,nrms_,&ierr_);
      // check that the rank deficiency was detected
      ASSERT_EQ(1, ierr_);
      // check norms
      _MT_ nrms_ref[_NV_];
      SUBR(mvec_norm2)(vec1_,nrms_ref,&ierr_);
      ASSERT_EQ(0,ierr_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_REAL_EQ(nrms_ref[i], nrms_[i]);
      }
      // check that we anyway got something orthogonal back
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      // the factor 2 in releps here is because otherwise fails the test by a fraction of releps
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.0*releps(vec1_));

      // check Q=V*B
      SUBR(mvec_times_sdMat)(-st::one(),vec1_,mat1_,st::one(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
      }
    }

  TEST_F(CLASSNAME, DISABLED_with_one_vectors) 
  {
    if (typeImplemented_)
    {
      SUBR(mvec_put_value)(vec1_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(svqb)(vec2_,mat1_,nrms_,&ierr_);
      // check that the rank deficiency was detected
      ASSERT_EQ(std::max(nvec_-1,0), ierr_);
      // check norms
      _MT_ nrms_ref[_NV_];
      SUBR(mvec_norm2)(vec1_,nrms_ref,&ierr_);
      ASSERT_EQ(0,ierr_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_REAL_EQ(nrms_ref[i], nrms_[i]);
      }
      // check that we anyway got something orthogonal back
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_DEB("B=\n");
      SUBR(sdMat_print)(mat1_,&ierr_);
#endif

      // check Q=V*B
      SUBR(mvec_times_sdMat)(-st::one(),vec1_,mat1_,st::one(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
    }
  }

