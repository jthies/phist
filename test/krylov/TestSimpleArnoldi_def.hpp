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
class CLASSNAME: public KernelTestWithSdMats<_ST_,_M_+BLOCK_SIZE1,_M_>,
                 public KernelTestWithVectors<_ST_,_N_,_M_+BLOCK_SIZE1>
{

  public:

    typedef KernelTestWithSdMats<_ST_,_M_+BLOCK_SIZE1,_M_> MTest;
    typedef KernelTestWithVectors<_ST_,_N_,_M_+BLOCK_SIZE1> VTest;

    //! mvec/sdMat sizes
    static const int n_=_N_;
    static const int m_=_M_;

    //! Set up routine.
    virtual void SetUp()
    {
      int iflag;
      MTest::SetUp();
      VTest::SetUp();

      if( typeImplemented_ )
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
        typeImplemented_ = typeImplemented_ && !globalTooSmall;
      }

      if( typeImplemented_ )
      {
        // read matrices
        SUBR(read_mat)("speye",nglob_,&Aeye_,&iflag_);
        ASSERT_EQ(0,iflag_);
#ifndef SKIP_ZERO_MAT
        SUBR(read_mat)("spzero",nglob_,&Azero_,&iflag_);
        ASSERT_EQ(0,iflag_);
#else
        Azero_=NULL;
#endif        
        SUBR(read_mat)("sprandn",nglob_,&Arand_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(read_mat)("sprandn_nodiag",nglob_,&Arand_nodiag_,&iflag_);
        ASSERT_EQ(0,iflag_);

        ASSERT_TRUE(Aeye_ != NULL);
#ifndef SKIP_ZERO_MAT        
        ASSERT_TRUE(Azero_ != NULL);
#endif
        ASSERT_TRUE(Arand_ != NULL);
        ASSERT_TRUE(Arand_nodiag_ != NULL);

        // wrap matrices in operators
        opAeye_ = new TYPE(op);
        ASSERT_TRUE(opAeye_ != NULL);
        opAzero_ = new TYPE(op);
        ASSERT_TRUE(opAzero_ != NULL);
        opArand_ = new TYPE(op);
        ASSERT_TRUE(opArand_ != NULL);
        opArand_nodiag_ = new TYPE(op);
        ASSERT_TRUE(opArand_nodiag_ != NULL);

        SUBR(op_wrap_sparseMat)(opAeye_,Aeye_,&iflag);
        ASSERT_EQ(0,iflag);
#ifndef SKIP_ZERO_MAT
        SUBR(op_wrap_sparseMat)(opAzero_,Azero_,&iflag);
        ASSERT_EQ(0,iflag);
#endif
        SUBR(op_wrap_sparseMat)(opArand_,Arand_,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(op_wrap_sparseMat)(opArand_nodiag_,Arand_nodiag_,&iflag);
        ASSERT_EQ(0,iflag);

        // setup views for needed vectors and sdMats
        v0_ = V_ = Vm_ = AV_ = AVm_ = VH_ = NULL;
        SUBR(mvec_view_block)(vec2_,&v0_,0,0,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_view_block)(vec1_,&V_,0,m_+BLOCK_SIZE1-1,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&iflag);
        ASSERT_EQ(0,iflag);
        stride_ = 1;
        SUBR(mvec_view_block)(vec1_,&Vm_,0,m_-1,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_view_block)(vec2_,&AV_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_view_block)(vec2_,&AVm_,0,m_-1,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_view_block)(vec3_,&VH_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&iflag);
        ASSERT_EQ(0,iflag);

        H_ = Hm_ = NULL;
        SUBR(sdMat_view_block)(mat1_,&H_,0,m_+BLOCK_SIZE1-1,0,m_-1,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(sdMat_view_block)(mat2_,&Hm_,0,m_-1,0,m_-1,&iflag);
        ASSERT_EQ(0,iflag);
      }
    }

    virtual void replaceMap(const_map_ptr_t map)
    {
      if( typeImplemented_ )
      {
        // delete vecs
        SUBR(mvec_delete)(v0_,&iflag_); v0_ = NULL;
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(V_,&iflag_); V_ = NULL;
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(Vm_,&iflag_); Vm_ = NULL;
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(AV_,&iflag_); AV_ = NULL;
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(AVm_,&iflag_); AVm_ = NULL;
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(VH_,&iflag_); VH_ = NULL;
        ASSERT_EQ(0,iflag_);

        VTest::replaceMap(map);

        // recreate vecs
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
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      int iflag;
      if( typeImplemented_ )
      {
        SUBR(mvec_delete)(v0_,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_delete)(V_, &iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_delete)(Vm_, &iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_delete)(AV_,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_delete)(AVm_,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvec_delete)(VH_,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(sdMat_delete)(H_,&iflag);
        ASSERT_EQ(0,iflag);
        SUBR(sdMat_delete)(Hm_,&iflag);
        ASSERT_EQ(0,iflag);

        delete opAeye_;
        delete opAzero_;
        delete opArand_;
        delete opArand_nodiag_;

        ASSERT_EQ(0,delete_mat(Aeye_));
        ASSERT_EQ(0,delete_mat(Azero_));
        ASSERT_EQ(0,delete_mat(Arand_));
        ASSERT_EQ(0,delete_mat(Arand_nodiag_));
      }

      VTest::TearDown();
      MTest::TearDown();
    }

    TYPE(op_ptr) opAeye_;
    TYPE(op_ptr) opAzero_;
    TYPE(op_ptr) opArand_;
    TYPE(op_ptr) opArand_nodiag_;

  protected:
    TYPE(sparseMat_ptr) Aeye_;
    TYPE(sparseMat_ptr) Azero_;
    TYPE(sparseMat_ptr) Arand_;
    TYPE(sparseMat_ptr) Arand_nodiag_;

    TYPE(mvec_ptr) v0_;
    TYPE(mvec_ptr) V_;
    _ST_ *V_vp_; lidx_t ldaV_, stride_;
    TYPE(mvec_ptr) Vm_;
    TYPE(mvec_ptr) AV_;
    TYPE(mvec_ptr) AVm_;
    TYPE(mvec_ptr) VH_;
    TYPE(sdMat_ptr) H_;
    TYPE(sdMat_ptr) Hm_;

    int delete_mat(TYPE(sparseMat_ptr) A)
    {
      if (A!=NULL)
        SUBR(sparseMat_delete)(A,&iflag_);
      return iflag_;
    }



    // ========================= the actual arnoldi test =========================
    void doArnoldiTest(TYPE(const_op_ptr) opA, int blockSize = 0)
    {
      int iflag;
      if( typeImplemented_ )
      {
        // run simple_arnoldi
        if( blockSize > 0 )
        {
          SUBR(simple_blockArnoldi)(opA,NULL,V_,NULL,NULL,H_,m_,blockSize,&iflag);
          ASSERT_EQ(0,iflag);
        }
        else
        {
          SUBR(simple_arnoldi)(opA,NULL,v0_,V_,NULL,NULL,H_,m_,&iflag);
          ASSERT_EQ(0,iflag);
        }

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&iflag);
        ASSERT_EQ(0,iflag);
        // calculate V(:,1:m+BLOCK_SIZE1)*H(1:m+BLOCK_SIZE1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&iflag);
        ASSERT_EQ(0,iflag);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&iflag);
        ASSERT_EQ(0,iflag);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&iflag);
        ASSERT_EQ(0,iflag);

        // check AV_' = AV_ - VH_ == 0
        for(int i = 0; i < _M_; i++)
          ASSERT_NEAR(mt::zero(),vnorm[i],100*mt::eps());
      }
    }

    // ========================= extended arnoldi test =========================
    // when one is interested in AV as well
    void doExtendedArnoldiTest(TYPE(const_op_ptr) opA, int blockSize = 0)
    {
      int iflag;
      if( typeImplemented_ )
      {
        // run simple_arnoldi
        if( blockSize > 0 )
        {
          SUBR(simple_blockArnoldi)(opA,NULL,V_,AVm_,NULL,H_,m_,blockSize,&iflag);
          ASSERT_EQ(0,iflag);
        }
        else
        {
          SUBR(simple_arnoldi)(opA,NULL,v0_,V_,AVm_,NULL,H_,m_,&iflag);
          ASSERT_EQ(0,iflag);
        }

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)200.*releps(V_));

        // check H = Vm'*AVm
        SUBR(sdMat_put_value)(mat2_, st::zero(), &iflag);
        ASSERT_EQ(0,iflag);
        SUBR(sdMat_get_block)(H_, Hm_, 0,m_-1, 0,m_-1, &iflag);
        ASSERT_EQ(0,iflag);
        SUBR(mvecT_times_mvec)(-st::one(),Vm_,AVm_,st::one(),Hm_, &iflag);
        ASSERT_EQ(0,iflag);
        ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,m_-1,m_-1,m_lda_,stride_,st::zero()), (MT)200.*releps(V_));

        // check A*Vm = AVm
        opA->apply(-st::one(),opA->A,Vm_,st::one(),AVm_,&iflag);
        ASSERT_EQ(0,iflag);
#ifdef PHIST_MVECS_ROW_MAJOR
        ASSERT_NEAR(mt::one(),ArrayEqual(vec2_vp_,m_+BLOCK_SIZE1,nloc_,lda_,stride_,st::zero()), (MT)200.*releps(V_));
#else
        ASSERT_NEAR(mt::one(),ArrayEqual(vec2_vp_,nloc_,m_+BLOCK_SIZE1,lda_,stride_,st::zero()), (MT)200.*releps(V_));
#endif

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&iflag);
        ASSERT_EQ(0,iflag);
        // calculate V(:,1:m+BLOCK_SIZE1)*H(1:m+BLOCK_SIZE1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&iflag);
        ASSERT_EQ(0,iflag);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&iflag);
        ASSERT_EQ(0,iflag);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&iflag);
        ASSERT_EQ(0,iflag);

        // check AV_' = AV_ - VH_ == 0
        for(int i = 0; i < _M_; i++)
          ASSERT_NEAR(mt::zero(),vnorm[i],100*mt::eps());
      }
    }
};


TEST_F(CLASSNAME, Aeye_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAeye_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doArnoldiTest(opAeye_, BLOCK_SIZE);
  }
}
#ifndef SKIP_ZERO_MAT
TEST_F(CLASSNAME, Azero_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAzero_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doArnoldiTest(opAzero_, BLOCK_SIZE);
  }
}
#endif
TEST_F(CLASSNAME, Arand_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doArnoldiTest(opArand_, BLOCK_SIZE);
  }
}

TEST_F(CLASSNAME, Arand_nodiag_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_nodiag_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doArnoldiTest(opArand_nodiag_, BLOCK_SIZE);
  }
}


TEST_F(CLASSNAME, extended_Aeye_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAeye_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doExtendedArnoldiTest(opAeye_, BLOCK_SIZE);
  }
}

#ifndef SKIP_ZERO_MAT
TEST_F(CLASSNAME, extended_Azero_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAzero_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doExtendedArnoldiTest(opAzero_, BLOCK_SIZE);
  }
}
#endif
TEST_F(CLASSNAME, extended_Arand_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doExtendedArnoldiTest(opArand_, BLOCK_SIZE);
  }
}

TEST_F(CLASSNAME, extended_Arand_nodiag_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_nodiag_->domain_map);
    ASSERT_EQ(0,iflag_);

    int iflag;
    SUBR(mvec_put_value)(v0_,st::one(),&iflag);
    ASSERT_EQ(0,iflag);
    doExtendedArnoldiTest(opArand_nodiag_, BLOCK_SIZE);
  }
}

