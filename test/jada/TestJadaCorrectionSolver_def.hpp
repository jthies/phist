#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>
{

  public:
    typedef KernelTestWithVectors<_ST_,_N_,_NV_> VTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NVP_> VProjTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;


    /*! Set up routine.
    */
    virtual void SetUp()
    {
      VTest::SetUp();
      MTest::SetUp();

      if (typeImplemented_)
      {
        SUBR(mvec_create)(&q_,map_,_NVP_,&ierr_);
        ASSERT_EQ(0,ierr_);
        sigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
          sigma_[i] = st::rand();

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&ierr_);
        ASSERT_EQ(0,ierr_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_QR)(q_,Rtmp,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(sdMat_delete)(Rtmp,&ierr_);
        ASSERT_EQ(0,ierr_);

        SUBR(read_mat)(MATNAME,nglob_,&A_,&ierr_);
        ASSERT_EQ(0,ierr_);
        ASSERT_TRUE(A_ != NULL);
        opA_ = new TYPE(op);
        SUBR(op_wrap_crsMat)(opA_, A_, &ierr_);
        ASSERT_EQ(0,ierr_);

        jdOp_ = new TYPE(op);
        SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,jdOp_,&ierr_);
        ASSERT_EQ(0,ierr_);

        // setup system to solve, exact x and A*x
        SUBR(mvec_random)(vec2_,&ierr_);
        ASSERT_EQ(0,ierr_);
        // x needs to be in q^orth
        SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_times_sdMat)(-st::one(),q_,mat1_,st::one(),vec2_,&ierr_);
        ASSERT_EQ(0,ierr_);
        jdOp_->apply(st::one(),jdOp_->A,vec2_,st::zero(),vec3_,&ierr_);
        ASSERT_EQ(0,ierr_);
      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      MTest::TearDown();
      VTest::TearDown();
      if (typeImplemented_)
      {
        SUBR(jadaOp_delete)(jdOp_,&ierr_);
        ASSERT_EQ(0,ierr_);
        delete jdOp_;
        delete opA_;
        SUBR(crsMat_delete)(A_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(q_,&ierr_);
        ASSERT_EQ(0,ierr_);
        delete[] sigma_;
      }
    }

    TYPE(crsMat_ptr) A_;
    TYPE(op_ptr) opA_;
    TYPE(op_ptr) jdOp_;
    TYPE(mvec_ptr) q_;
    _ST_* sigma_;
};


  TEST_F(CLASSNAME, create_and_delete)
  {
    if( typeImplemented_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      SUBR(jadaCorrectionSolver_create)(solver, 3, map_, 10, &ierr_);
      ASSERT_EQ(0, ierr_);
      SUBR(jadaCorrectionSolver_delete)(solver, &ierr_);
      ASSERT_EQ(0, ierr_);
    }
  }

