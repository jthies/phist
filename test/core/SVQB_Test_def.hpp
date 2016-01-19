#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,2>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_> 
  {

public:

  typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,2> VTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_-1> VTest_m_minus_1;
  typedef KernelTestWithVectors<_ST_,_N_,1> VTest_1;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_> MTest;

  static void SetUpTestCase()
  {
    VTest::SetUpTestCase();
    MTest::SetUpTestCase();
  }

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    VTest::SetUp();
    MTest::SetUp();
    if (typeImplemented_ && !problemTooSmall_)
      {
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
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
      SUBR(mvec_to_device)(vec2_,&iflag_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    VTest::TearDown();
    MTest::TearDown();
    }

  static void TearDownTestCase()
  {
    MTest::TearDownTestCase();
    VTest::TearDownTestCase();
  }



  _MT_ nrms_[_NV_];
};

  TEST_F(CLASSNAME, with_random_vectors) 
    {
    if (typeImplemented_ && !problemTooSmall_)
      {
     //PrintVector(*cout,"SVQB_Test V",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      SUBR(svqb)(vec2_,mat1_,nrms_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // doing a SVQB decomp must not relocate data:
      ASSERT_EQ(true,MTest::pointerUnchanged(mat1_,mat1_vp_,m_lda_));
      ASSERT_EQ(true,VTest::pointerUnchanged(vec2_,vec2_vp_,lda_));
      //PrintVector(*cout,"SVQB_Test Q",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      /*
      PHIST_SOUT(PHIST_DEBUG,"B matrix:\n");
      SUBR(sdMat_print)(mat1_,&iflag_);
      */
      // check norms
      _MT_ nrms_ref[_NV_];
      SUBR(mvec_norm2)(vec1_,nrms_ref,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(nrms_ref[i], nrms_[i],100*nrms_ref[i]*mt::eps());
      }
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));

      // check Q=V*B
      SUBR(mvec_times_sdMat)(-st::one(),vec1_,mat1_,st::one(),vec2_,&iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
      }
    }

  TEST_F(CLASSNAME, with_rank_deficiency) 
    {
    if (typeImplemented_ && !problemTooSmall_)
      {
      if (nvec_==1)
        {
        SUBR(mvec_put_value)(vec1_,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);
        }
      else
        {
        SUBR(mvec_random)(vec1_,&iflag_);
        SUBR(mvec_from_device)(vec1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // set last two columns to same vector
        for (int i=0;i<stride_*nloc_;i+=stride_)
          {
#ifdef PHIST_MVECS_ROW_MAJOR
          vec1_vp_[(nvec_-1)+i*lda_] = vec1_vp_[(nvec_-2)+i*lda_];
#else
          vec1_vp_[(nvec_-1)*lda_+i] = vec1_vp_[(nvec_-2)*lda_+i];
#endif
          }
        SUBR(mvec_to_device)(vec1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        }
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(svqb)(vec2_,mat1_,nrms_,&iflag_);
      // check that the rank deficiency was detected
      ASSERT_EQ(1, iflag_);
      // check norms
      _MT_ nrms_ref[_NV_];
      iflag_ = 0;
      SUBR(mvec_norm2)(vec1_,nrms_ref,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(nrms_ref[i], nrms_[i],100*nrms_ref[i]*mt::eps());
      }
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // SVQB will return a block with the first <dim0> columns zero,
      // and the remaining orthonormal
      _ST_* vec2_col1=NULL;
#ifdef PHIST_MVECS_ROW_MAJOR
      vec2_col1 = vec2_vp_+1;
#else
      vec2_col1 = vec2_vp_+lda_;
#endif

      // the factor 2 in releps here is because otherwise fails the test by a fraction of releps
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,1,lda_,stride_,st::zero(),vflag_), 100*mt::eps());
      if (nvec_>1)
      {
        ASSERT_NEAR(mt::one(),VTest_m_minus_1::ColsAreNormalized(vec2_col1,nloc_,lda_,stride_,mpi_comm_),(MT)100.0*releps(vec1_));
        ASSERT_NEAR(mt::one(),VTest_m_minus_1::ColsAreOrthogonal(vec2_col1,nloc_,lda_,stride_,mpi_comm_),(MT)100.0*releps(vec1_));
      }

      // check Q=V*B
      SUBR(mvec_times_sdMat)(-st::one(),vec1_,mat1_,st::one(),vec2_,&iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
      }
    }

  TEST_F(CLASSNAME, with_one_vectors) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      SUBR(mvec_put_value)(vec1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(svqb)(vec2_,mat1_,nrms_,&iflag_);
      // check that the rank deficiency was detected
      ASSERT_EQ(std::max(nvec_-1,0), iflag_);
      // check norms
      _MT_ nrms_ref[_NV_];
      SUBR(mvec_norm2)(vec1_,nrms_ref,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(nrms_ref[i], nrms_[i],10*nrms_ref[i]*mt::eps());
      }
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // SVQB will return a block with the first <dim0> columns zero,
      // and the remaining orthonormal
      _ST_* vec2_col1=NULL;
#ifdef PHIST_MVECS_ROW_MAJOR
      vec2_col1 = vec2_vp_+nvec_-1;
#else
      vec2_col1 = vec2_vp_+(nvec_-1)*lda_;
#endif
      // the factor 2 in releps here is because otherwise fails the test by a fraction of releps
      ASSERT_EQ(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_-1,lda_,stride_,st::zero(),vflag_));
      if (nvec_>1)
      {
        ASSERT_NEAR(mt::one(),VTest_1::ColsAreNormalized(vec2_col1,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      }

      // check Q=V*B
      SUBR(mvec_times_sdMat)(-st::one(),vec1_,mat1_,st::one(),vec2_,&iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
    }
  }

