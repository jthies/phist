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
        SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,mat1_,vec1_,NULL,NULL,jdOp_,&ierr_);
        ASSERT_EQ(0,ierr_);

        // setup system to solve, exact x and A*x
        SUBR(mvec_random)(vec2_,&ierr_);
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
      TYPE(jadaInnerGmresState_ptr) state[_NV_];
      SUBR(jadaInnerGmresStates_create)(state, _NV_, map_, _MAXBAS_, &ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(jadaInnerGmresStates_delete)(state, _NV_, &ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, reset_and_updateSol_single_without_iteration)
  {
    if( typeImplemented_ )
    {
      TYPE(jadaInnerGmresState_ptr) state[_NV_];
      SUBR(jadaInnerGmresStates_create)(state, _NV_, map_, _MAXBAS_, &ierr_);
      ASSERT_EQ(0,ierr_);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);

        SUBR(jadaInnerGmresState_reset)(state[i], y_i, x_i, &ierr_);
        ASSERT_EQ(0,ierr_);
      }

      // don't let it iterate here!

      // check the result (we have given the solution as initial guess!)
      TYPE(mvec_ptr) Ax_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec1_,&Ax_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);

        _MT_ resNorm;
        SUBR(jadaInnerGmresStates_updateSol)(&state[i], 1, x_i, Ax_i, &resNorm, &ierr_);
        ASSERT_EQ(0,ierr_);

        // can't know the residual norm, yet
        ASSERT_REAL_EQ(-1,resNorm);
      }
      // delete views
      SUBR(mvec_delete)(x_i,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(Ax_i,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(y_i,&ierr_);
      ASSERT_EQ(0,ierr_);

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      // check Ax = A*x
      SUBR(crsMat_times_mvec)(-st::one(),A_,vec2_,st::one(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif


      SUBR(jadaInnerGmresStates_delete)(state, _NV_, &ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }
