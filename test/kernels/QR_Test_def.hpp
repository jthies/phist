#ifndef CLASSNAME
#error "file not included correctly"
#endif

#if !defined(MVECS_VIEWED)|| !defined(SDMATS_VIEWED)
#error "file not included correctly"
#endif

#ifdef TEST_MVEC_QR
#define MVEC_QR SUBR(mvec_QR)
#elif defined(TEST_CHOL_QR)
#define MVEC_QR SUBR(chol_QR)
#endif

/*! Test the kernel function mvec_QR. As this is an 'optional' kernel, all tests pass if they encounter
    a return value of PHIST_NOT_IMPLEMENTED in mvec_QR. The orthog routine (core lib) still works if 
    mvec_QR is not available by using our own implementation of CholQR. 
 
 The same source file is used to generate the tests for chol_QR, which has the same signature and 
 semantics as mvec_QR.
 */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_,MVECS_VIEWED,2>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_,SDMATS_VIEWED> 
{

public:

  typedef KernelTestWithVectors<_ST_,_N_,_NV_,MVECS_VIEWED,2> VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_,SDMATS_VIEWED> MTest;

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

    if(typeImplemented_ && !problemTooSmall_)
    {

      // disable the test because TSQR will not work.
      // This situation occurs if we have a small matrix (_N_=25, say)
      // and many Q vectors (e.g. 10) with multiple MPI procs.
      int globalTooSmall = _N_ < _NV_;
#ifdef PHIST_HAVE_MPI
      int localTooSmall = nloc_ < _NV_;
      iflag_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
      ASSERT_EQ(0,iflag_);
#endif
      problemTooSmall_ = globalTooSmall != 0;
    }

    if (typeImplemented_ && !problemTooSmall_)
      {
      SUBR(mvec_random)(vec1_,&iflag_);
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
      ASSERT_EQ(0,iflag_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    MTest::TearDown();
    VTest::TearDown();
    }

  static void TearDownTestCase()
  {
    MTest::TearDownTestCase();
    VTest::TearDownTestCase();
  }

};

// this is in fact an MvecTest, but here we have the useful 'ColsAreNormalized' test
TEST_F(CLASSNAME,mvec_normalize)
{
  if (typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _MT_ norms[nvec_];
    for (int j=0;j<nvec_;j++) norms[j]=mt::one()/(_MT_)(j+1);
    SUBR(mvec_normalize)(vec1_,norms,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),ColsAreNormalized(vec1_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec2_));

    for (int i=0;i<nloc_;i++)
      for (int j=0;j<nvec_;j++)
      {
        vec1_vp_[VIDX(i,j,lda_)]*=norms[j];
      }
    SUBR(mvec_to_device)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_REAL_EQ(mt::one(), MvecsEqual(vec1_,vec2_));
  }
}

  TEST_F(CLASSNAME, with_random_vectors) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
//      PrintVector(*cout,"QR_Test V",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      MVEC_QR(vec2_,mat1_,&iflag_);
#ifdef TEST_MVEC_QR
      if (iflag_==PHIST_NOT_IMPLEMENTED) return;
#endif
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // doing a QR decomp must not relocate data:
      ASSERT_EQ(true,MTest::pointerUnchanged(mat1_,mat1_vp_,m_lda_));
      ASSERT_EQ(true,VTest::pointerUnchanged(vec2_,vec2_vp_,lda_));
//      PrintVector(*cout,"QR_Test Q",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));

      // check Q*R=V
      SUBR(mvec_times_sdMat)(-st::one(),vec2_,mat1_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
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
        ASSERT_EQ(0,iflag_);
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
      }
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      MVEC_QR(vec2_,mat1_,&iflag_);
#ifdef TEST_MVEC_QR
      if (iflag_==PHIST_NOT_IMPLEMENTED) return;
#endif
      // check that the rank deficiency was detected
      ASSERT_EQ(1, iflag_);
      iflag_ = 0;
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check that we anyway got something orthogonal back
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      // the factor 2 in releps here is because otherwise fails the test by a fraction of releps
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.0*releps(vec1_));

#if PHIST_OUTLEV>=PHIST_DEBUG
# if _N_<100
        PHIST_SOUT(PHIST_DEBUG,"original V:\n");
        SUBR(mvec_print)(vec1_,&iflag_);
        PHIST_SOUT(PHIST_DEBUG,"computed Q:\n");
        SUBR(mvec_print)(vec2_,&iflag_);
        PHIST_SOUT(PHIST_DEBUG,"computed R:\n");
        SUBR(sdMat_print)(mat1_,&iflag_);
# endif
#endif
      // check Q*R=V
      SUBR(mvec_times_sdMat)(-st::one(),vec2_,mat1_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
# if _N_<100
        PHIST_SOUT(PHIST_DEBUG,"explicit Q*R:\n");
        SUBR(mvec_print)(vec1_,&iflag_);
# endif
#endif
      ASSERT_NEAR(mt::one(), ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
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

      MVEC_QR(vec2_,mat1_,&iflag_);
#ifdef TEST_MVEC_QR
      if (iflag_==PHIST_NOT_IMPLEMENTED) return;
#endif
      // check that the rank deficiency was detected
      ASSERT_EQ(std::max(nvec_-1,0), iflag_);
      iflag_ = 0;
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check that we anyway got something orthogonal back
      ASSERT_NEAR(mt::one(),ColsAreNormalized(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
      ASSERT_NEAR(mt::one(),ColsAreOrthogonal(vec2_vp_,nloc_,lda_,stride_,mpi_comm_),(MT)100.*releps(vec1_));
#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_DEB("R=\n");
      SUBR(sdMat_print)(mat1_,&iflag_);
#endif
      // check Q*R=V
      SUBR(mvec_times_sdMat)(-st::one(),vec2_,mat1_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(), ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_),sqrt(mt::eps()));
    }
  }

