/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#if !defined(CLASSNAME) || !defined(INNER_SOLVTYPE)
#error "file not included correctly."
#endif

/*! Test fixture. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME>,
                 public virtual KernelTestWithMassMat<_ST_,_N_>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_M_,0,3>,
                 public virtual KernelTestWithSdMats<_ST_,_M_,_M_>,
                 public JadaTestWithOpts
{

  public:
  
    typedef KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME> ATest;
    typedef KernelTestWithMassMat<_ST_,_N_>               BTest;
    typedef KernelTestWithVectors<_ST_,_N_,_M_,0,3>       VTest;
    typedef KernelTestWithSdMats<_ST_,_M_,_M_>            MTest;
    typedef TestWithType<_MT_>                            RTest;

    static constexpr _ST_ BmatScaling_=(_ST_)8;
    static const int  blockSize_=_K_;
    static const int  nEig_=_M_-_K_+1;
    static const int minBas_=2*_M_;
    static const int maxBas_=minBas_+8*blockSize_;
    static const int maxIters_=20*(int)((maxBas_-minBas_+1)/blockSize_);
    static constexpr _MT_ convTol_=std::sqrt(mt::eps());

    static void SetUpTestCase()
    {
      int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
      ATest::SetUpTestCase(sparseMatCreateFlag);
      BTest::SetUpTestCase(context_,phist::testing::PHIST_TG_PREFIX(idfunc),(void*)&BmatScaling_);
      VTest::SetUpTestCase();
      MTest::SetUpTestCase();
    }

    /*! Set up routine.
    */
    virtual void SetUp()
    {
      ATest::SetUp();
      BTest::SetUp();
      VTest::SetUp();
      MTest::SetUp();
      JadaTestWithOpts::SetUp();

      jadaOpts_.innerSolvType=INNER_SOLVTYPE;
      jadaOpts_.minBas=minBas_;
      jadaOpts_.maxBas=maxBas_;
      jadaOpts_.convTol=convTol_;
      jadaOpts_.maxIters=maxIters_;
      jadaOpts_.blockSize=blockSize_;

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
        SUBR(linearOp_wrap_sparseMat)(opA_, A_, &iflag_);
        ASSERT_EQ(0,iflag_);

        opB_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opB_, B_, &iflag_);
        ASSERT_EQ(0,iflag_);

        opAB_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat_pair)(opAB_, A_, B_, &iflag_);
        ASSERT_EQ(0,iflag_);


      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      if (typeImplemented_ && !problemTooSmall_)
      {
        if( opA_ != NULL )
          delete opA_;
        opA_ = NULL;

        if( opB_ != NULL )
          delete opB_;
        opB_ = NULL;

        if( opAB_ != NULL )
          delete opAB_;
        opAB_ = NULL;
      }
      
      MTest::TearDown();
      VTest::TearDown();
      BTest::TearDown();
      ATest::TearDown();
    }

    static void TearDownTestCase()
    {
      MTest::TearDownTestCase();
      VTest::TearDownTestCase();
      BTest::TearDownTestCase();
      ATest::TearDownTestCase();
    }
/*
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
*/
    TYPE(linearOp_ptr) opA_ = NULL, opB_=NULL, opAB_=NULL;
};


     constexpr _ST_ CLASSNAME::BmatScaling_;
     const int  CLASSNAME::blockSize_;
     const int  CLASSNAME::nEig_;
     const int CLASSNAME::minBas_;
     const int CLASSNAME::maxBas_;
     const int CLASSNAME::maxIters_;
     constexpr _MT_ CLASSNAME::convTol_;



  TEST_F(CLASSNAME, with_and_without_B)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(mvec_ptr)  Q = VTest::vec1_;
      TYPE(sdMat_ptr) R = MTest::mat1_;
      int nEig=nEig_;
      int nIter;
      _CT_ ev[_M_];
      _MT_ resNorm[_M_];
      SUBR(subspacejada)(opA_, NULL, jadaOpts_,
                     Q, R, ev, resNorm, &nEig, &nIter, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nEig,nEig_);
      
      // run the same thing with B=i*alpha and make sure the eigenvalues are just scaled by alpha
      TYPE(mvec_ptr) Qb=VTest::vec2_;
      TYPE(sdMat_ptr) Rb=MTest::mat2_;
      _CT_ evb[nEig_+blockSize_-1];
      SUBR(subspacejada)(opAB_, opB_, jadaOpts_,
                     Qb, Rb, evb, resNorm, &nEig, &nIter, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nEig,nEig_);
      
      for (int i=0; i<nEig; i++) evb[i]/BmatScaling_;
      // check real parts
      ASSERT_NEAR(mt::one(),RTest::ArraysEqual((_MT_*)ev,(_MT_*)evb,nEig,1,1,2),jadaOpts_.convTol);
      // check imaginary parts
      ASSERT_NEAR(mt::one(),RTest::ArraysEqual(((_MT_*)ev)+1,((_MT_*)evb)+1,nEig,1,1,2),(_MT_)(12*jadaOpts_.convTol));
    
    }
  }
  
