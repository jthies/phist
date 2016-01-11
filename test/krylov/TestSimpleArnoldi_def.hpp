#ifndef CLASSNAME
#error "file not included correctly"
#endif

#ifdef BLOCK_SIZE1
#undef BLOCK_SIZE1
#endif
#if BLOCK_SIZE > 0
#define BLOCK_SIZE1 BLOCK_SIZE
#else
#define BLOCK_SIZE1 1
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithSdMats<_ST_,_M_+BLOCK_SIZE1,_M_>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_M_+BLOCK_SIZE1,0,3>
{

  public:
    typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
    typedef KernelTestWithVectors<_ST_,_N_,_M_+BLOCK_SIZE1,0,3> VTest;
    typedef KernelTestWithSdMats<_ST_,_M_+BLOCK_SIZE1,_M_> MTest;

    //! mvec/sdMat sizes
    static const int n_=_N_;
    static const int m_=_M_;

    static void SetUpTestCase()
    {
      SparseMatTest::SetUpTestCase();
      VTest::SetUpTestCase();
      MTest::SetUpTestCase();
    }

    //! Set up routine.
    virtual void SetUp()
    {
      SparseMatTest::SetUp();
      MTest::SetUp();
      VTest::SetUp();

      if( typeImplemented_ && !problemTooSmall_ )
      {

        // disable the test because TSQR will not work.
        // This situation occurs if we have a small matrix (_N_=25, say)
        // and many Q vectors (e.g. 10) with multiple MPI procs.
        int globalTooSmall = _N_ < BLOCK_SIZE1;
#ifdef PHIST_HAVE_MPI
        int localTooSmall = nloc_ < BLOCK_SIZE1;
        iflag_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        ASSERT_EQ(0,iflag_);
#endif
        problemTooSmall_ = globalTooSmall != 0;
      }

      if( typeImplemented_ && !problemTooSmall_ )
      {
        // wrap matrices in operators
        opA_ = new TYPE(op);
        ASSERT_TRUE(opA_ != NULL);

        SUBR(op_wrap_sparseMat)(opA_,A_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // setup views for needed vectors and sdMats
        v0_ = V_ = Vm_ = AV_ = AVm_ = VH_ = NULL;
        SUBR(mvec_view_block)(vec2_,&v0_,0,0,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec1_,&V_,0,m_+BLOCK_SIZE1-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&iflag_);
        ASSERT_EQ(0,iflag_);
        stride_ = 1;
        SUBR(mvec_view_block)(vec1_,&Vm_,0,m_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec2_,&AV_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec2_,&AVm_,0,m_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&VH_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&iflag_);
        ASSERT_EQ(0,iflag_);

        H_ = Hm_ = NULL;
        SUBR(sdMat_view_block)(mat1_,&H_,0,m_+BLOCK_SIZE1-1,0,m_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_view_block)(mat2_,&Hm_,0,m_-1,0,m_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      if( typeImplemented_ && !problemTooSmall_ )
      {
        SUBR(mvec_delete)(v0_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(V_, &iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(Vm_, &iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(AV_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(AVm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(VH_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_delete)(H_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_delete)(Hm_,&iflag_);
        ASSERT_EQ(0,iflag_);

        delete opA_;
        opA_ = NULL;
      }

      MTest::TearDown();
      VTest::TearDown();
      SparseMatTest::TearDown();
    }

    static void TearDownTestCase()
    {
      MTest::TearDownTestCase();
      VTest::TearDownTestCase();
      SparseMatTest::TearDownTestCase();
    }

    TYPE(op_ptr) opA_;

  protected:
    TYPE(mvec_ptr) v0_;
    TYPE(mvec_ptr) V_;
    _ST_ *V_vp_; lidx_t ldaV_, stride_;
    TYPE(mvec_ptr) Vm_;
    TYPE(mvec_ptr) AV_;
    TYPE(mvec_ptr) AVm_;
    TYPE(mvec_ptr) VH_;
    TYPE(sdMat_ptr) H_;
    TYPE(sdMat_ptr) Hm_;


    // ========================= the actual arnoldi test =========================
    void doArnoldiTest(TYPE(const_op_ptr) opA, int blockSize = 0)
    {
      if( typeImplemented_ && !problemTooSmall_ )
      {
        // run simple_arnoldi
        if( blockSize > 0 )
        {
          SUBR(simple_blockArnoldi)(opA,NULL,V_,NULL,NULL,H_,m_,blockSize,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        else
        {
          SUBR(simple_arnoldi)(opA,NULL,v0_,V_,NULL,NULL,H_,m_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // calculate V(:,1:m+BLOCK_SIZE1)*H(1:m+BLOCK_SIZE1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&iflag_);
        ASSERT_EQ(0,iflag_);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&iflag_);
        ASSERT_EQ(0,iflag_);

        // check AV_' = AV_ - VH_ == 0
        for(int i = 0; i < _M_; i++)
        {
          ASSERT_NEAR(mt::zero(),vnorm[i],(MT)200.*releps(V_));
        }
      }
    }

    // ========================= extended arnoldi test =========================
    // when one is interested in AV as well
    void doExtendedArnoldiTest(TYPE(const_op_ptr) opA, int blockSize = 0)
    {
      if( typeImplemented_ && !problemTooSmall_ )
      {
        // run simple_arnoldi
        if( blockSize > 0 )
        {
          SUBR(simple_blockArnoldi)(opA,NULL,V_,AVm_,NULL,H_,m_,blockSize,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        else
        {
          SUBR(simple_arnoldi)(opA,NULL,v0_,V_,AVm_,NULL,H_,m_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));

        // check H = Vm'*AVm
        SUBR(sdMat_put_value)(mat2_, st::zero(), &iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_get_block)(H_, Hm_, 0,m_-1, 0,m_-1, &iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvecT_times_mvec)(-st::one(),Vm_,AVm_,st::one(),Hm_, &iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,m_-1,m_-1,m_lda_,stride_,st::zero()), (MT)200.*releps(V_));

        // check A*Vm = AVm
        opA->apply(-st::one(),opA->A,Vm_,st::one(),AVm_,&iflag_);
        ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
        ASSERT_NEAR(mt::one(),ArrayEqual(vec2_vp_,m_+BLOCK_SIZE1,nloc_,lda_,stride_,st::zero()), (MT)200.*releps(V_));
#else
        ASSERT_NEAR(mt::one(),ArrayEqual(vec2_vp_,nloc_,m_+BLOCK_SIZE1,lda_,stride_,st::zero()), (MT)200.*releps(V_));
#endif

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // calculate V(:,1:m+BLOCK_SIZE1)*H(1:m+BLOCK_SIZE1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&iflag_);
        ASSERT_EQ(0,iflag_);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&iflag_);
        ASSERT_EQ(0,iflag_);

        // check AV_' = AV_ - VH_ == 0
        for(int i = 0; i < _M_; i++)
          ASSERT_NEAR(mt::zero(),vnorm[i],(MT)200.*releps(V_));
      }
    }
};


#if MATNAME == MATNAME_speye
TEST_F(CLASSNAME, Aeye_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

#if MATNAME == MATNAME_spzero
TEST_F(CLASSNAME, Azero_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

#if MATNAME == MATNAME_sprandn
TEST_F(CLASSNAME, Arand_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

#if MATNAME == MATNAME_sprandn_nodiag
TEST_F(CLASSNAME, Arand_nodiag_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif


#if MATNAME == MATNAME_speye
TEST_F(CLASSNAME, extended_Aeye_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doExtendedArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

#if MATNAME == MATNAME_spzero
TEST_F(CLASSNAME, extended_Azero_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doExtendedArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

#if MATNAME == MATNAME_sprandn
TEST_F(CLASSNAME, extended_Arand_v0ones) 
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doExtendedArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

#if MATNAME == MATNAME_sprandn_nodiag
TEST_F(CLASSNAME, extended_Arand_nodiag_v0ones)
{
  if( typeImplemented_ && !problemTooSmall_)
  {
    SUBR(mvec_put_value)(v0_,st::one(),&iflag_);
    ASSERT_EQ(0,iflag_);
    doExtendedArnoldiTest(opA_, BLOCK_SIZE);
  }
}
#endif

