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

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif
      return;

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


  TEST_F(CLASSNAME, reset_and_updateSol_vector_without_iteration)
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
      // delete views
      SUBR(mvec_delete)(x_i,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(y_i,&ierr_);
      ASSERT_EQ(0,ierr_);

      // don't let it iterate here!

      // check the result (we have given the solution as initial guess!)
      _MT_ resNorm[_NV_];
      SUBR(jadaInnerGmresStates_updateSol)(state, _NV_, vec2_, vec3_, resNorm, &ierr_);
      ASSERT_EQ(0,ierr_);

      // can't know the residual norm, yet
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_REAL_EQ(-1,resNorm[i]);
      }

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


  TEST_F(CLASSNAME, iterate_with_exact_initial_guess)
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
      // delete views
      SUBR(mvec_delete)(x_i,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(y_i,&ierr_);
      ASSERT_EQ(0,ierr_);

      // call iterate
      int nIter = 0;
      SUBR(jadaInnerGmresStates_iterate)(jdOp_,state, _NV_, &nIter, &ierr_);
      ASSERT_EQ(0,ierr_);
      // only one iteration should be needed!
      ASSERT_EQ(1,nIter);
      // all systems should be marked as converged
      for(int i = 0; i < _NV_; i++)
      {
        // all systems did 1 iteration
        ASSERT_EQ(1,state[i]->totalIter);
        // all systems converged
        ASSERT_EQ(0,state[i]->ierr);
      }

      SUBR(crsMat_times_mvec)(st::one(),A_,vec2_,st::zero(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // check the result (we have given the solution as initial guess!)
      _MT_ resNorm[_NV_];
      SUBR(jadaInnerGmresStates_updateSol)(state, _NV_, vec2_, vec1_, resNorm, &ierr_);
      ASSERT_EQ(0,ierr_);

      // can't know the residual norm, yet
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(0,resNorm[i],10*VTest::releps());
      }

      // check Ax = A*x
      SUBR(crsMat_times_mvec)(-st::one(),A_,vec2_,st::one(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif


      SUBR(jadaInnerGmresStates_delete)(state, _NV_, &ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, iterate_with_single_exact_initial_guess)
  {
    if( typeImplemented_ )
    {
      TYPE(jadaInnerGmresState_ptr) state[_NV_];
      SUBR(jadaInnerGmresStates_create)(state, _NV_, map_, _MAXBAS_, &ierr_);
      ASSERT_EQ(0,ierr_);

      int exactGuessAt = std::min(1,_NV_-1);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&ierr_);
        ASSERT_EQ(0,ierr_);

        if( i != exactGuessAt )
        {
          SUBR(mvec_put_value)(x_i,st::zero(),&ierr_);
          ASSERT_EQ(0,ierr_);
        }
        SUBR(jadaInnerGmresState_reset)(state[i], y_i, x_i, &ierr_);
        ASSERT_EQ(0,ierr_);
      }
      // delete views
      SUBR(mvec_delete)(x_i,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(y_i,&ierr_);
      ASSERT_EQ(0,ierr_);

      // call iterate
      int nIter = 0;
      SUBR(jadaInnerGmresStates_iterate)(jdOp_,state, _NV_, &nIter, &ierr_);
      ASSERT_EQ(0,ierr_);
      // only one iteration should be needed!
      ASSERT_EQ(1,nIter);
      // at least one system converged
      ASSERT_EQ(0,state[exactGuessAt]->ierr);
      for(int i = 0; i < _NV_; i++)
      {
        // all systems did 1 iteration
        ASSERT_EQ(1,state[i]->totalIter);
        // other systems may have converged
        ASSERT_TRUE(state[i]->ierr == 0 || state[i]->ierr == 1);
      }

      SUBR(crsMat_times_mvec)(st::one(),A_,vec2_,st::zero(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // check the result
      _MT_ resNorm[_NV_];
      SUBR(jadaInnerGmresStates_updateSol)(state, _NV_, vec2_, vec1_, resNorm, &ierr_);
      ASSERT_EQ(0,ierr_);

      // resnorm should be set now
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_OUT(PHIST_INFO,"resNorm[%d] = %8.4e\n", resNorm[i]);
        ASSERT_TRUE(resNorm[i] >= mt::zero());
      }
      ASSERT_NEAR(mt::zero(), resNorm[exactGuessAt], 10*VTest::releps());

      // check Ax = A*x
      SUBR(crsMat_times_mvec)(-st::one(),A_,vec2_,st::one(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif


      // now check the result: vec3[exactGuessAt] = jdOp_(vec2)[exactGuessAt]
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(&vec3_vp_[exactGuessAt],1,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(&vec3_vp_[exactGuessAt*lda_],nloc_,1,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(jadaInnerGmresStates_delete)(state, _NV_, &ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, iterate_a_bit)
  {
    if( typeImplemented_ )
    {
      TYPE(jadaInnerGmresState_ptr) state[_NV_];
      SUBR(jadaInnerGmresStates_create)(state, _NV_, map_, _MAXBAS_, &ierr_);
      ASSERT_EQ(0,ierr_);

      // setup system with y = 0.3*jdOp^2(x) - 0.9*jdOp(x)
      {
        TYPE(mvec_ptr) tmp;
        SUBR(mvec_create)(&tmp,map_,_NV_,&ierr_);
        ASSERT_EQ(0,ierr_);

        SUBR(mvec_add_mvec)(st::one(),vec3_,st::zero(),tmp,&ierr_);
        ASSERT_EQ(0,ierr_);
        jdOp_->apply((_ST_)0.3*st::one(),jdOp_->A,tmp,(_ST_)-0.9*st::one(),vec3_,&ierr_);
        ASSERT_EQ(0,ierr_);

        SUBR(mvec_delete)(tmp,&ierr_);
        ASSERT_EQ(0,ierr_);
      }
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

        state[i]->tol = 100*VTest::releps();
      }
      // delete views
      SUBR(mvec_delete)(x_i,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(y_i,&ierr_);
      ASSERT_EQ(0,ierr_);

      // call iterate
      int nIter = 0;
      SUBR(jadaInnerGmresStates_iterate)(jdOp_,state, _NV_, &nIter, &ierr_);
      ASSERT_TRUE(ierr_ == 0 || ierr_ == 1);
      for(int i = 0; i < _NV_; i++)
      {
        // all systems did nIter iteration
        ASSERT_EQ(nIter,state[i]->totalIter);
      }

      SUBR(crsMat_times_mvec)(st::one(),A_,vec2_,st::zero(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // check the result (we have given the solution as initial guess!)
      _MT_ resNorm[_NV_];
      SUBR(jadaInnerGmresStates_updateSol)(state, _NV_, vec2_, vec1_, resNorm, &ierr_);
      ASSERT_EQ(0,ierr_);

      // check Ax = A*x
      SUBR(crsMat_times_mvec)(-st::one(),A_,vec2_,st::one(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      // now check the result: vec3 = jdOp_(vec2)
      //jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&ierr_);
      //ASSERT_EQ(0,ierr_);
//#ifdef PHIST_KERNEL_LIB_FORTRAN
      //ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
//#else
      //ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
//#endif


      SUBR(jadaInnerGmresStates_delete)(state, _NV_, &ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


