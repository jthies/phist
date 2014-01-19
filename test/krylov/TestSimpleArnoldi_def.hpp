#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_M_+1,_M_>,
                 public KernelTestWithVectors<_ST_,_N_,_M_+1>
{

  public:

    typedef KernelTestWithSdMats<_ST_,_M_+1,_M_> MTest;
    typedef KernelTestWithVectors<_ST_,_N_,_M_+1> VTest;

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
        // read matrices
        SUBR(read_mat)("speye",nglob_,&Aeye_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(read_mat)("spzero",nglob_,&Azero_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(read_mat)("sprandn",nglob_,&Arand_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(read_mat)("sprandn_nodiag",nglob_,&Arand_nodiag_,&ierr_);
        ASSERT_EQ(0,ierr_);

        ASSERT_TRUE(Aeye_ != NULL);
        ASSERT_TRUE(Azero_ != NULL);
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
        SUBR(op_wrap_crsMat)(opAzero_,Azero_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(op_wrap_crsMat)(opArand_,Arand_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(op_wrap_crsMat)(opArand_nodiag_,Arand_nodiag_,&ierr);
        ASSERT_EQ(0,ierr);

        // setup views for needed vectors and sdMats
        v0_ = V_ = Vm_ = AV_ = VH_ = NULL;
        SUBR(mvec_view_block)(vec2_,&v0_,0,0,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec1_,&V_,0,m_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&ierr);
        ASSERT_EQ(0,ierr);
        stride_ = 1;
        SUBR(mvec_view_block)(vec1_,&Vm_,0,m_-1,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec2_,&AV_,1,m_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(mvec_view_block)(vec3_,&VH_,1,m_,&ierr);
        ASSERT_EQ(0,ierr);

        H_ = NULL;
        SUBR(sdMat_view_block)(mat1_,&H_,0,m_,0,m_-1,&ierr);
        ASSERT_EQ(0,ierr);
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
        SUBR(mvec_delete)(VH_,&ierr);
        ASSERT_EQ(0,ierr);
        SUBR(sdMat_delete)(H_,&ierr);
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
    _ST_ *V_vp_; int ldaV_, stride_;
    TYPE(mvec_ptr) Vm_;
    TYPE(mvec_ptr) AV_;
    TYPE(mvec_ptr) VH_;
    TYPE(sdMat_ptr) H_;

    int delete_mat(TYPE(crsMat_ptr) A)
    {
      if (A!=NULL)
        SUBR(crsMat_delete)(A,&ierr_);
      return ierr_;
    }



    // ========================= the actual arnoldi test =========================
    void doArnoldiTest(TYPE(const_op_ptr) opA)
    {
      int ierr;
      if( typeImplemented_ )
      {
        // run simple_arnoldi
        SUBR(simple_arnoldi)(opA,NULL,v0_,V_,NULL,NULL,H_,m_,&ierr);
        ASSERT_EQ(0,ierr);

        // check orthogonality of V_
        ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)50.*releps(V_));
        ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_),(MT)50.*releps(V_));

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        // calculate V(:,1:m+1)*H(1:m+1,1:m)
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
    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opAeye_);
  }
}

TEST_F(CLASSNAME, Azero_v0ones) 
{
  if( typeImplemented_)
  {
    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opAzero_);
  }
}

TEST_F(CLASSNAME, Arand_v0ones) 
{
  if( typeImplemented_)
  {
    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opArand_);
  }
}

TEST_F(CLASSNAME, Arand_nodiag_v0ones) 
{
  if( typeImplemented_)
  {
    int ierr;
    SUBR(mvec_put_value)(v0_,st::one(),&ierr);
    ASSERT_EQ(0,ierr);
    doArnoldiTest(opArand_nodiag_);
  }
}

