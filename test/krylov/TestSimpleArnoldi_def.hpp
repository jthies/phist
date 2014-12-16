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
      int ierr;
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
        ierr_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        ASSERT_EQ(0,ierr_);
#endif
        typeImplemented_ = typeImplemented_ && !globalTooSmall;
      }

      if( typeImplemented_ )
      {
        // read matrices
        SUBR(read_mat)("speye",nglob_,&Aeye_,&ierr_);
        ASSERT_EQ(0,ierr_);
#ifndef SKIP_ZERO_MAT
        SUBR(read_mat)("spzero",nglob_,&Azero_,&ierr_);
        ASSERT_EQ(0,ierr_);
#else
        Azero_=NULL;
#endif        
        SUBR(read_mat)("sprandn",nglob_,&Arand_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(read_mat)("sprandn_nodiag",nglob_,&Arand_nodiag_,&ierr_);
        ASSERT_EQ(0,ierr_);

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

        SUBR(op_wrap_crsMat)(opAeye_,Aeye_,&ierr);
        ASSERT_EQ(0,ierr);
#ifndef SKIP_ZERO_MAT
        SUBR(op_wrap_crsMat)(opAzero_,Azero_,&ierr);
        ASSERT_EQ(0,ierr);
#endif
        SUBR(op_wrap_crsMat)(opArand_,Arand_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(op_wrap_crsMat)(opArand_nodiag_,Arand_nodiag_,&ierr);
        ASSERT_EQ(0,ierr);

        // setup views for needed vectors and sdMats
        v0_ = V_ = Vm_ = AV_ = AVm_ = VH_ = NULL;
        SUBR(mvec_view_block)(vec2_,&v0_,0,0,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec1_,&V_,0,m_+BLOCK_SIZE1-1,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&ierr);
        ASSERT_EQ(0,ierr);
        stride_ = 1;
        SUBR(mvec_view_block)(vec1_,&Vm_,0,m_-1,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec2_,&AV_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec2_,&AVm_,0,m_-1,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec3_,&VH_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&ierr);
        ASSERT_EQ(0,ierr);

        H_ = Hm_ = NULL;
        SUBR(sdMat_view_block)(mat1_,&H_,0,m_+BLOCK_SIZE1-1,0,m_-1,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(sdMat_view_block)(mat2_,&Hm_,0,m_-1,0,m_-1,&ierr);
        ASSERT_EQ(0,ierr);
      }
    }

    virtual void replaceMap(const_map_ptr_t map)
    {
      if( typeImplemented_ )
      {
        // delete vecs
        SUBR(mvec_delete)(v0_,&ierr_); v0_ = NULL;
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(V_,&ierr_); V_ = NULL;
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(Vm_,&ierr_); Vm_ = NULL;
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(AV_,&ierr_); AV_ = NULL;
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(AVm_,&ierr_); AVm_ = NULL;
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(VH_,&ierr_); VH_ = NULL;
        ASSERT_EQ(0,ierr_);

        VTest::replaceMap(map);

        // recreate vecs
        SUBR(mvec_view_block)(vec2_,&v0_,0,0,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec1_,&V_,0,m_+BLOCK_SIZE1-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&ierr_);
        ASSERT_EQ(0,ierr_);
        stride_ = 1;
        SUBR(mvec_view_block)(vec1_,&Vm_,0,m_-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec2_,&AV_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec2_,&AVm_,0,m_-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec3_,&VH_,BLOCK_SIZE1,m_+BLOCK_SIZE1-1,&ierr_);
        ASSERT_EQ(0,ierr_);
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      int ierr;
      if( typeImplemented_ )
      {
        SUBR(mvec_delete)(v0_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_delete)(V_, &ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_delete)(Vm_, &ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_delete)(AV_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_delete)(AVm_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_delete)(VH_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(sdMat_delete)(H_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(sdMat_delete)(Hm_,&ierr);
        ASSERT_EQ(0,ierr);

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
    TYPE(crsMat_ptr) Aeye_;
    TYPE(crsMat_ptr) Azero_;
    TYPE(crsMat_ptr) Arand_;
    TYPE(crsMat_ptr) Arand_nodiag_;

    TYPE(mvec_ptr) v0_;
    TYPE(mvec_ptr) V_;
    _ST_ *V_vp_; lidx_t ldaV_, stride_;
    TYPE(mvec_ptr) Vm_;
    TYPE(mvec_ptr) AV_;
    TYPE(mvec_ptr) AVm_;
    TYPE(mvec_ptr) VH_;
    TYPE(sdMat_ptr) H_;
    TYPE(sdMat_ptr) Hm_;

    int delete_mat(TYPE(crsMat_ptr) A)
    {
      if (A!=NULL)
        SUBR(crsMat_delete)(A,&ierr_);
      return ierr_;
    }



    // ========================= the actual arnoldi test =========================
    void doArnoldiTest(TYPE(const_op_ptr) opA, int blockSize = 0)
    {
      int ierr;
      if( typeImplemented_ )
      {
        // run simple_arnoldi
        if( blockSize > 0 )
        {
          SUBR(simple_blockArnoldi)(opA,NULL,V_,NULL,NULL,H_,m_,blockSize,&ierr);
          ASSERT_EQ(0,ierr);
        }
        else
        {
          SUBR(simple_arnoldi)(opA,NULL,v0_,V_,NULL,NULL,H_,m_,&ierr);
          ASSERT_EQ(0,ierr);
        }

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)50.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)50.*releps(V_));

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        // calculate V(:,1:m+BLOCK_SIZE1)*H(1:m+BLOCK_SIZE1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&ierr);
        ASSERT_EQ(0,ierr);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&ierr);
        ASSERT_EQ(0,ierr);

        // check AV_' = AV_ - VH_ == 0
        for(int i = 0; i < _M_; i++)
          ASSERT_NEAR(mt::zero(),vnorm[i],100*mt::eps());
      }
    }

    // ========================= extended arnoldi test =========================
    // when one is interested in AV as well
    void doExtendedArnoldiTest(TYPE(const_op_ptr) opA, int blockSize = 0)
    {
      int ierr;
      if( typeImplemented_ )
      {
        // run simple_arnoldi
        if( blockSize > 0 )
        {
          SUBR(simple_blockArnoldi)(opA,NULL,V_,AVm_,NULL,H_,m_,blockSize,&ierr);
          ASSERT_EQ(0,ierr);
        }
        else
        {
          SUBR(simple_arnoldi)(opA,NULL,v0_,V_,AVm_,NULL,H_,m_,&ierr);
          ASSERT_EQ(0,ierr);
        }

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)50.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)50.*releps(V_));

        // check H = Vm'*AVm
        SUBR(sdMat_put_value)(mat2_, st::zero(), &ierr);
        ASSERT_EQ(0,ierr);
        SUBR(sdMat_get_block)(H_, Hm_, 0,m_-1, 0,m_-1, &ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvecT_times_mvec)(-st::one(),Vm_,AVm_,st::one(),Hm_, &ierr);
        ASSERT_EQ(0,ierr);
        ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,m_-1,m_-1,m_lda_,stride_,st::zero()), (MT)50.*releps(V_));

        // check A*Vm = AVm
        opA->apply(-st::one(),opA->A,Vm_,st::one(),AVm_,&ierr);
        ASSERT_EQ(0,ierr);
#ifdef PHIST_MVECS_ROW_MAJOR
        ASSERT_NEAR(mt::one(),ArrayEqual(vec2_vp_,m_+BLOCK_SIZE1,nloc_,lda_,stride_,st::zero()), (MT)50.*releps(V_));
#else
        ASSERT_NEAR(mt::one(),ArrayEqual(vec2_vp_,nloc_,m_+BLOCK_SIZE1,lda_,stride_,st::zero()), (MT)50.*releps(V_));
#endif

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        // calculate V(:,1:m+BLOCK_SIZE1)*H(1:m+BLOCK_SIZE1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&ierr);
        ASSERT_EQ(0,ierr);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&ierr);
        ASSERT_EQ(0,ierr);

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
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opAeye_, BLOCK_SIZE);
  }
}
#ifndef SKIP_ZERO_MAT
TEST_F(CLASSNAME, Azero_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAzero_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opAzero_, BLOCK_SIZE);
  }
}
#endif
TEST_F(CLASSNAME, Arand_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opArand_, BLOCK_SIZE);
  }
}

TEST_F(CLASSNAME, Arand_nodiag_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_nodiag_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opArand_nodiag_, BLOCK_SIZE);
  }
}


TEST_F(CLASSNAME, extended_Aeye_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAeye_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doExtendedArnoldiTest(opAeye_, BLOCK_SIZE);
  }
}

#ifndef SKIP_ZERO_MAT
TEST_F(CLASSNAME, extended_Azero_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opAzero_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doExtendedArnoldiTest(opAzero_, BLOCK_SIZE);
  }
}
#endif
TEST_F(CLASSNAME, extended_Arand_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doExtendedArnoldiTest(opArand_, BLOCK_SIZE);
  }
}

TEST_F(CLASSNAME, extended_Arand_nodiag_v0ones) 
{
  if( typeImplemented_)
  {
    replaceMap(opArand_nodiag_->domain_map);
    ASSERT_EQ(0,ierr_);

    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doExtendedArnoldiTest(opArand_nodiag_, BLOCK_SIZE);
  }
}

