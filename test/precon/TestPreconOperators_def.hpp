/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME>,
                 public virtual KernelTestWithSparseMat<_ST_,_N_,_N_,PRECNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,4>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>,
                 public JadaTestWithOpts
{

  public:
    typedef KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME> ATest;
    typedef KernelTestWithSparseMat<_ST_,_N_,_N_,PRECNAME> PTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,4> VTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;

    static void SetUpTestCase()
    {
      int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
      ATest::SetUpTestCase(sparseMatCreateFlag);
      PTest::SetUpTestCase(sparseMatCreateFlag);
      VTest::SetUpTestCase();
      MTest::SetUpTestCase();
    }

    /*! Set up routine.
    */
    virtual void SetUp()
    {
      ATest::SetUp();
      PTest::SetUp();
      VTest::SetUp();
      MTest::SetUp();
      JadaTestWithOpts::SetUp();

      jadaOpts_.innerSolvType=phist_GMRES;
      jadaOpts_.maxBas=_MAXBAS_;

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
        opA_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opA_, ATest::A_, &iflag_);
        ASSERT_EQ(0,iflag_);

        opP_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opP_, PTest::A_, &iflag_);
        ASSERT_EQ(0,iflag_);
        
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

        jdOp_ = new TYPE(linearOp);
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
      if (typeImplemented_ && !problemTooSmall_)
      {
        SUBR(jadaOp_delete)(jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( jdOp_ != NULL )
          delete jdOp_;
        jdOp_ = NULL;
        if( opA_ != NULL )
          delete opA_;
        opA_ = NULL;
        if( opP_ != NULL )
          delete opP_;
        opP_ = NULL;
        SUBR(mvec_delete)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( sigma_ != NULL )
          delete[] sigma_;
        sigma_ = NULL;
        if( negSigma_ != NULL )
          delete[] negSigma_;
        negSigma_ = NULL;
      }
      MTest::TearDown();
      VTest::TearDown();
      ATest::TearDown();
      PTest::TearDown();
    }

    static void TearDownTestCase()
    {
      MTest::TearDownTestCase();
      VTest::TearDownTestCase();
      ATest::TearDownTestCase();
      PTest::TearDownTestCase();
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

    TYPE(linearOp_ptr) opA_ = NULL, opP_ = NULL;
    TYPE(linearOp_ptr) jdOp_ = NULL, jdPrec_ = NULL;
    TYPE(mvec_ptr) q_ = NULL;
    _ST_* sigma_ = NULL;
    _ST_* negSigma_ = NULL;
};

  TEST_F(CLASSNAME, wrap_Ainv_as_USER_PRECON)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
        // now we have our approximate inverse as a linearOp, we can wrap it once more to get
        // a valid phist preconditioner. This object adds an 'update' function etc, which should
        // return with an error. It also takes care of checking for apply_shifted and falling back
        // to apply if it is NULL.
        TYPE(linearOp) userPrec;
        // we explicitly disable the 'apply_shifted' function of our approximate inverse
        opP_->apply_shifted=NULL;
        // note that we do not need to pass in any matrix because we don't provide an 'update' function
        SUBR(precon_create)(&userPrec,NULL,sigma_[0],NULL,NULL,NULL,
                "user_defined",NULL,opP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // aplying the preconditioner should be the same as applying the original Ainv matrix
        _ST_ alpha=st::prand(),beta=st::prand();
        SUBR(mvec_random)(vec1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_random)(vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sparseMat_times_mvec)(alpha,PTest::A_,vec1_,beta,vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(precon_apply)(alpha,PTest::A_,vec1_,beta,vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_),10*VTest::releps());
    }
  }

  TEST_F(CLASSNAME, apply_jadaPrec)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      // apply wrapped projected preconditioner and check orthogonality of result
    }
  }


  TEST_F(CLASSNAME, precon_improves_krylov_solver)
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
        SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, &sigma_[i],res_i, NULL, &tol[i], 200, t_i, 1, 0, 0, &iflag_);
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


  TEST_F(CLASSNAME, jadaCorrectionSolver)
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
      SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, sigma_,vec3_, NULL, tol, 200, vec2_, 1, 0, 0, &iflag_);
      ASSERT_EQ(0, iflag_);

      // check all solutions
      checkResiduals(tol);

      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }

  

  TEST_F(CLASSNAME, subspacejada)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
    }
  }
