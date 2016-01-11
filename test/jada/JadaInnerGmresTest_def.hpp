#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>
{

  public:
    typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NV_> VTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;

    static void SetUpTestCase()
    {
      SparseMatTest::SetUpTestCase();
      MTest::SetUpTestCase();
      VTest::SetUpTestCase();
    }

    /*! Set up routine.
    */
    virtual void SetUp()
    {
      SparseMatTest::SetUp();
      VTest::SetUp();
      MTest::SetUp();

      if( typeImplemented_ && !problemTooSmall_ )
      {
        // disable the test because TSQR will not work.
        // This situation occurs if we have a small matrix (_N_=25, say)
        // and many Q vectors (e.g. 10) with multiple MPI procs.
        int globalTooSmall = _N_ <= std::min(_NVP_,_NV_);
#ifdef PHIST_HAVE_MPI
        int localTooSmall = nloc_ <= std::min(_NVP_,_NV_);
        iflag_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        ASSERT_EQ(0,iflag_);
#endif
        problemTooSmall_ = globalTooSmall != 0;
      }

      if (typeImplemented_ && !problemTooSmall_)
      {
        opA_ = new TYPE(op);
        SUBR(op_wrap_sparseMat)(opA_, A_, &iflag_);
        ASSERT_EQ(0,iflag_);

        PHISTTEST_MVEC_CREATE(&q_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        sigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
          sigma_[i] = st::prand();

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        int rankQ=0;
        SUBR(orthog)(NULL,q_,NULL,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(sdMat_delete)(Rtmp,&iflag_);
        ASSERT_EQ(0,iflag_);

        jdOp_ = new TYPE(op);
        SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // setup system to solve, exact x and A*x
        SUBR(mvec_random)(vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // x needs to be in q^orth
        SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_times_sdMat)(-st::one(),q_,mat1_,st::one(),vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        jdOp_->apply(st::one(),jdOp_->A,vec2_,st::zero(),vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      if (typeImplemented_ && !problemTooSmall_)
      {
        SUBR(jadaOp_delete)(jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( jdOp_ != NULL )
          delete jdOp_;
        if( opA_ != NULL )
          delete opA_;
        SUBR(mvec_delete)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( sigma_ != NULL )
          delete[] sigma_;
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


    TYPE(op_ptr) opA_ = NULL;
    TYPE(op_ptr) jdOp_ = NULL;
    TYPE(mvec_ptr) q_ = NULL;
    _ST_* sigma_ = NULL;
};


  TEST_F(CLASSNAME, create_and_delete)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(blockedGMRESstate_ptr) state[_NV_];
      SUBR(blockedGMRESstates_create)(state, _NV_, map_, _MAXBAS_, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(blockedGMRESstates_delete)(state, _NV_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, reset_and_updateSol_single_without_iteration)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(blockedGMRESstate_ptr) state[_NV_];
      SUBR(blockedGMRESstates_create)(state, _NV_, map_, _MAXBAS_, &iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(blockedGMRESstate_reset)(state[i], y_i, x_i, &iflag_);
        ASSERT_EQ(0,iflag_);
      }

      // don't let it iterate here!
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);

        _MT_ resNorm;
        SUBR(blockedGMRESstates_updateSol)(&state[i], 1, x_i, &resNorm, false, &iflag_);
        ASSERT_EQ(0,iflag_);

        // can't know the residual norm, yet
        ASSERT_REAL_EQ(-1,resNorm);
      }
      // delete views
      SUBR(mvec_delete)(x_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(y_i,&iflag_);
      ASSERT_EQ(0,iflag_);

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(blockedGMRESstates_delete)(state, _NV_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, reset_and_updateSol_vector_without_iteration)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(blockedGMRESstate_ptr) state[_NV_];
      SUBR(blockedGMRESstates_create)(state, _NV_, map_, _MAXBAS_, &iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(blockedGMRESstate_reset)(state[i], y_i, x_i, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
      // delete views
      SUBR(mvec_delete)(x_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(y_i,&iflag_);
      ASSERT_EQ(0,iflag_);

      // don't let it iterate here!

      // check the result (we have given the solution as initial guess!)
      _MT_ resNorm[_NV_];
      SUBR(blockedGMRESstates_updateSol)(state, _NV_, vec2_, resNorm, false, &iflag_);
      ASSERT_EQ(0,iflag_);

      // can't know the residual norm, yet
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_REAL_EQ(-1,resNorm[i]);
      }

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(blockedGMRESstates_delete)(state, _NV_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, iterate_with_exact_initial_guess)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(blockedGMRESstate_ptr) state[_NV_];
      SUBR(blockedGMRESstates_create)(state, _NV_, map_, _MAXBAS_, &iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(blockedGMRESstate_reset)(state[i], y_i, x_i, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
      // delete views
      SUBR(mvec_delete)(x_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(y_i,&iflag_);
      ASSERT_EQ(0,iflag_);

      // call iterate
      int nIter = 0;
      SUBR(blockedGMRESstates_iterate)(jdOp_,state, _NV_, &nIter, true, &iflag_);
      ASSERT_EQ(0,iflag_);
      // only one iteration should be needed!
      ASSERT_EQ(1,nIter);
      // all systems should be marked as converged
      for(int i = 0; i < _NV_; i++)
      {
        // all systems did 1 iteration
        ASSERT_EQ(1,state[i]->totalIter);
        // all systems converged
        ASSERT_EQ(0,state[i]->status);
      }

 
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec2_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check the result (we have given the solution as initial guess!)
      _MT_ resNorm[_NV_];
      SUBR(blockedGMRESstates_updateSol)(state, _NV_, vec2_, resNorm, false, &iflag_);
      ASSERT_EQ(0,iflag_);

      // residual didn't change, so the relative residual should be one
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(1,resNorm[i],10*VTest::releps());
      }

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif


      SUBR(blockedGMRESstates_delete)(state, _NV_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, iterate_with_single_exact_initial_guess)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(blockedGMRESstate_ptr) state[_NV_];
      SUBR(blockedGMRESstates_create)(state, _NV_, map_, _MAXBAS_, &iflag_);
      ASSERT_EQ(0,iflag_);

      int exactGuessAt = std::min(1,_NV_-1);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);

        if( i != exactGuessAt )
        {
          SUBR(mvec_put_value)(x_i,st::zero(),&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        SUBR(blockedGMRESstate_reset)(state[i], y_i, x_i, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
      // delete views
      SUBR(mvec_delete)(x_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(y_i,&iflag_);
      ASSERT_EQ(0,iflag_);

      // call iterate
      int nIter = 0;
      SUBR(blockedGMRESstates_iterate)(jdOp_,state, _NV_, &nIter, true, &iflag_);
      ASSERT_EQ(0,iflag_);
      // only one iteration should be needed!
      ASSERT_EQ(1,nIter);
      // at least one system converged
      ASSERT_EQ(0,state[exactGuessAt]->status);
      for(int i = 0; i < _NV_; i++)
      {
        // all systems did 1 iteration
        ASSERT_EQ(1,state[i]->totalIter);
        // other systems may have converged
        ASSERT_TRUE(state[i]->status == 0 || state[i]->status == 1);
      }

      SUBR(sparseMat_times_mvec)(st::one(),A_,vec2_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check the result
      _MT_ resNorm[_NV_];
      SUBR(blockedGMRESstates_updateSol)(state, _NV_, vec2_, resNorm, false, &iflag_);
      ASSERT_EQ(0,iflag_);

      // resnorm should be set now
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO,"resNorm[%d] = %8.4e\n", i, resNorm[i]);
        ASSERT_TRUE(resNorm[i] >= mt::zero());
      }
      ASSERT_NEAR(mt::one(), resNorm[exactGuessAt], 10*VTest::releps());

      // now check the result: vec3[exactGuessAt] = jdOp_(vec2)[exactGuessAt]
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(&vec3_vp_[exactGuessAt],1,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(&vec3_vp_[exactGuessAt*lda_],nloc_,1,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(blockedGMRESstates_delete)(state, _NV_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, iterate_a_bit)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(blockedGMRESstate_ptr) state[_NV_];
      SUBR(blockedGMRESstates_create)(state, _NV_, map_, _MAXBAS_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // setup system with y = 0.3*jdOp^2(x) - 0.9*jdOp(x)
      {
        TYPE(mvec_ptr) tmp;
        PHISTTEST_MVEC_CREATE(&tmp,map_,_NV_,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(mvec_add_mvec)(st::one(),vec3_,st::zero(),tmp,&iflag_);
        ASSERT_EQ(0,iflag_);
        jdOp_->apply((_ST_)0.3*st::one(),jdOp_->A,tmp,(_ST_)-0.9*st::one(),vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(mvec_delete)(tmp,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      SUBR(mvec_put_value)(vec2_, st::one(), &iflag_);

      TYPE(mvec_ptr) x_i = NULL;
      TYPE(mvec_ptr) y_i = NULL;
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_,&x_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(vec3_,&y_i,i,i,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(blockedGMRESstate_reset)(state[i], y_i, x_i, &iflag_);
        ASSERT_EQ(0,iflag_);

        state[i]->tol = 100*VTest::releps();
      }
      // delete views
      SUBR(mvec_delete)(x_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(y_i,&iflag_);
      ASSERT_EQ(0,iflag_);

      // calculate initial residual norm
      SUBR(mvec_add_mvec)(st::one(),vec3_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // calculate the residual norm
      _MT_ initialResNorm[_NV_];
      SUBR(mvec_norm2)(vec1_,initialResNorm,&iflag_);

      // call iterate
      int nIter = 0;
      SUBR(blockedGMRESstates_iterate)(jdOp_,state, _NV_, &nIter, true, &iflag_);
      ASSERT_TRUE(iflag_ == 0 || iflag_ == 1);
      for(int i = 0; i < _NV_; i++)
      {
        // all systems did nIter iteration
        ASSERT_EQ(nIter,state[i]->totalIter);
      }

      SUBR(sparseMat_times_mvec)(st::one(),A_,vec2_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check the result (we have given the solution as initial guess!)
      _MT_ resNorm[_NV_];
      SUBR(blockedGMRESstates_updateSol)(state, _NV_, vec2_, resNorm, false, &iflag_);
      ASSERT_EQ(0,iflag_);

      // now check the result: vec3 = jdOp_(vec2)
      jdOp_->apply(-st::one(),jdOp_->A,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // calculate the residual norm
      _MT_ explicitResNorm[_NV_];
      SUBR(mvec_norm2)(vec3_,explicitResNorm,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHIST_SOUT(PHIST_INFO,"est. rel. residual norms:");
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO,"\t%8.4e", resNorm[i]);
      }
      PHIST_SOUT(PHIST_INFO,"\nexp. rel. residual norms:");
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO,"\t%8.4e", explicitResNorm[i]/initialResNorm[i]);
      }
      PHIST_SOUT(PHIST_INFO,"\n");
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_LT(resNorm[i], mt::one());
        ASSERT_LT(explicitResNorm[i], initialResNorm[i]);
        ASSERT_NEAR(explicitResNorm[i]/initialResNorm[i], resNorm[i], 100*VTest::releps());
      }
//#ifdef PHIST_MVECS_ROW_MAJOR
      //ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
//#else
      //ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
//#endif


      SUBR(blockedGMRESstates_delete)(state, _NV_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

