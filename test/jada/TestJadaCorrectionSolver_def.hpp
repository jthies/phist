#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>,
                 public JadaTestWithOpts
{

  public:
    typedef KernelTestWithVectors<_ST_,_N_,_NV_> VTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;


    /*! Set up routine.
    */
    virtual void SetUp()
    {
      VTest::SetUp();
      MTest::SetUp();
      JadaTestWithOpts::SetUp();

      jadaOpts_.innerSolvType=GMRES;
      jadaOpts_.maxBas=_MAXBAS_;

      if (typeImplemented_ && !problemTooSmall_)
      {
        SUBR(read_mat)(MATNAME,comm_,nglob_,&A_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_TRUE(A_ != NULL);
        opA_ = new TYPE(op);
        SUBR(op_wrap_sparseMat)(opA_, A_, &iflag_);
        ASSERT_EQ(0,iflag_);

        replaceMap(opA_->domain_map);

        PHISTTEST_MVEC_CREATE(&q_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        sigma_ = new _ST_[_NV_];
        negSigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
        {
          // there are hopefully no eigenvalues in this region so the matrix doesn't get nearly singular
          sigma_[i] = (_ST_)30*st::one() + (_ST_)5*st::prand();
          negSigma_[i] = -sigma_[i];
        }

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        int rankQ = 0;
        SUBR(orthog)(NULL,q_,NULL,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(sdMat_delete)(Rtmp,&iflag_);
        ASSERT_EQ(0,iflag_);

        jdOp_ = new TYPE(op);
        SUBR(jadaOp_create)(opA_,NULL,q_,NULL,negSigma_,_NV_,jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // setup system to solve, exact x and A*x
        SUBR(mvec_random)(vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // x needs to be in q^orth
        SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_times_sdMat)(-st::one(),q_,mat1_,st::one(),vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        _MT_ tmp[_NV_];
        SUBR(mvec_normalize)(vec2_, tmp, &iflag_);
        ASSERT_EQ(0,iflag_);
        jdOp_->apply(st::one(),jdOp_->A,vec2_,st::zero(),vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // as we are only interested in the direction of vec3_, scale it to one
        SUBR(mvec_normalize)(vec3_, tmp, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      MTest::TearDown();
      VTest::TearDown();
      if (typeImplemented_ && !problemTooSmall_)
      {
        SUBR(jadaOp_delete)(jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( jdOp_ != NULL )
          delete jdOp_;
        if( opA_ != NULL )
          delete opA_;
        SUBR(sparseMat_delete)(A_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( sigma_ != NULL )
          delete[] sigma_;
        if( negSigma_ != NULL )
          delete[] negSigma_;
      }
    }

    void checkResiduals(_MT_ tol[_NV_])
    {
      // the calculated approximations should be stored in vec2, the JD-residual in vec3_

      // the solution should be scaled to one
      _MT_ solutionNorm[_NV_];
      SUBR(mvec_norm2)(vec2_, solutionNorm, &iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(mt::one(), solutionNorm[i], 100*VTest::releps());
      }

      // we cannot directly compare the vectors because the solution is scaled
      jdOp_->apply(st::one(),jdOp_->A,vec2_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // rather check that vec1_ and vec3_ point in the same direction
      _MT_ tmp[_NV_];
      SUBR(mvec_normalize)(vec1_, tmp, &iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ dot[_NV_];
      SUBR(mvec_dot_mvec)(vec3_, vec1_, dot, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_vadd_mvec)(dot, vec3_, -st::one(), vec1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ resNorm[_NV_];
      SUBR(mvec_norm2)(vec1_, resNorm, &iflag_);
      ASSERT_EQ(0,iflag_);
      // check that the resnorm is actually near the required norm
      PHIST_SOUT(PHIST_INFO, "required tol.:");
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO, "\t%8.4e", tol[i]);
      }
      PHIST_SOUT(PHIST_INFO, "\nachieved tol.:");
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO, "\t%8.4e", resNorm[i]);
      }
      PHIST_SOUT(PHIST_INFO, "\n");
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_LT(resNorm[i], 20*tol[i]);
      }
    }

    TYPE(sparseMat_ptr) A_ = NULL;
    TYPE(op_ptr) opA_ = NULL;
    TYPE(op_ptr) jdOp_ = NULL;
    TYPE(mvec_ptr) q_ = NULL;
    _ST_* sigma_ = NULL;
    _ST_* negSigma_ = NULL;
};


  TEST_F(CLASSNAME, create_and_delete)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      jadaOpts_.blockSize=3;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);
      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }

  TEST_F(CLASSNAME, selftest)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      std::vector<_MT_> tol(_NV_, VTest::releps());
      checkResiduals(&tol[0]);
    }
  }

  TEST_F(CLASSNAME, single)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      jadaOpts_.blockSize=1;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);

      TYPE(mvec_ptr) t_i = NULL;
      TYPE(mvec_ptr) res_i = NULL;
      _MT_ tol[_NV_];
      for(int i = 0; i < _NV_; i++)
      {
        SUBR(mvec_view_block)(vec2_, &t_i, i, i, &iflag_);
        ASSERT_EQ(0, iflag_);
        SUBR(mvec_view_block)(vec3_, &res_i, i, i, &iflag_);
        ASSERT_EQ(0, iflag_);

        // create some random tolerance
        tol[i] = exp((_MT_)-8*mt::one() + (_MT_)4*mt::prand());

        SUBR(mvec_put_value)(t_i, st::zero(), &iflag_);
        ASSERT_EQ(0, iflag_);

        // run
        SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, &sigma_[i],res_i, NULL, &tol[i], 200, t_i, 1, 0, &iflag_);
        ASSERT_EQ(0, iflag_);
      }

      // check all solutions
      checkResiduals(tol);

      SUBR(mvec_delete)(res_i, &iflag_);
      ASSERT_EQ(0, iflag_);
      SUBR(mvec_delete)(t_i, &iflag_);
      ASSERT_EQ(0, iflag_);
      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }


  TEST_F(CLASSNAME, all_at_once)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      jadaOpts_.blockSize=_NV_;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);

      _MT_ tol[_NV_];
      for(int i = 0; i < _NV_; i++)
      {
        // create some random tolerance
        tol[i] = (_MT_)exp((_MT_)-8*mt::one() + (_MT_)4*mt::prand());
      }

      SUBR(mvec_put_value)(vec2_, st::zero(), &iflag_);
      ASSERT_EQ(0, iflag_);

      // run
      SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, sigma_,vec3_, NULL, tol, 200, vec2_, 1, 0, &iflag_);
      ASSERT_EQ(0, iflag_);

      // check all solutions
      checkResiduals(tol);

      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }


  TEST_F(CLASSNAME, one_after_another)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      jadaOpts_.blockSize=1;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);

      _MT_ tol[_NV_];
      for(int i = 0; i < _NV_; i++)
      {
        // create some random tolerance
        tol[i] = exp((_MT_)-8*mt::one() + (_MT_)4*mt::prand());
      }

      SUBR(mvec_put_value)(vec2_, st::zero(), &iflag_);
      ASSERT_EQ(0, iflag_);

      // run
      SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, sigma_, vec3_, NULL, tol, 200, vec2_, 1, 0, &iflag_);
      ASSERT_EQ(0, iflag_);

      // check all solutions
      checkResiduals(tol);

      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }


  TEST_F(CLASSNAME, pipelined_2)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      jadaOpts_.blockSize=2;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);

      _MT_ tol[_NV_];
      for(int i = 0; i < _NV_; i++)
      {
        // create some random tolerance
        tol[i] = (_MT_)exp((_MT_)-8*mt::one() + (_MT_)4*mt::prand());
      }

      SUBR(mvec_put_value)(vec2_, st::zero(), &iflag_);
      ASSERT_EQ(0, iflag_);

      // run
      SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, sigma_, vec3_, NULL, tol, 200, vec2_, 1, 0, &iflag_);
      ASSERT_EQ(0, iflag_);

      // check all solutions
      checkResiduals(tol);

      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }


  TEST_F(CLASSNAME, pipelined_4)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      jadaOpts_.blockSize=4;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);

      _MT_ tol[_NV_];
      for(int i = 0; i < _NV_; i++)
      {
        // create some random tolerance
        tol[i] = (_MT_)exp((_MT_)-8*mt::one() + (_MT_)4*mt::prand());
      }

      SUBR(mvec_put_value)(vec2_, st::zero(), &iflag_);
      ASSERT_EQ(0, iflag_);

      // run
      SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, sigma_, vec3_, NULL, tol, 200, vec2_, 1, 0, &iflag_);
      ASSERT_EQ(0, iflag_);

      // check all solutions
      checkResiduals(tol);

      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }


