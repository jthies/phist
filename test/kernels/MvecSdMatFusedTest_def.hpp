/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly"
#endif
#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif

#ifdef USE_VIEWS
#undef USE_VIEWS
#endif
#if _USE_VIEWS_V_ || _USE_VIEWS_W_ || _USE_VIEWS_M_ || _USE_VIEWS_N_
#define USE_VIEWS 1
#else
#define USE_VIEWS 0
#endif

/*! Test fixure. */
class CLASSNAME: public virtual TestWithType< _ST_ >,
                 public virtual KernelTestWithMap<_N_>
{

public:

  class VTest : public KernelTestWithVectors<_ST_,_N_,_M_,_USE_VIEWS_V_,2,1> {
    public: virtual void TestBody(){}
  };
  class WTest : public KernelTestWithVectors<_ST_,_N_,_K_,_USE_VIEWS_W_,2,2> {
    public: virtual void TestBody(){}
  };
  class MTest : public KernelTestWithSdMats<_ST_,_M_,_K_,_USE_VIEWS_M_,1> {
    public: virtual void TestBody(){}
  };
  class NTest : public KernelTestWithSdMats<_ST_,_K_,_K_,_USE_VIEWS_N_,2> {
    public: virtual void TestBody(){}
  };

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  static const int k_=_K_;
  
  //! V is n x m
  TYPE(mvec_ptr) V1_ = NULL, V2_ = NULL;

  //! W is n x k
  TYPE(mvec_ptr) W1_ = NULL, W2_ = NULL;

  //! M is m x k
  TYPE(sdMat_ptr) M1_ = NULL, M2_ = NULL;

  //! N is k x k
  TYPE(sdMat_ptr) N1_ = NULL, N2_ = NULL;
  
  _ST_ *V1_vp_,*V2_vp_,*W1_vp_,*W2_vp_,*M1_vp_,*M2_vp_,*N1_vp_,*N2_vp_;
  
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  phist_lidx ldaV1_,ldaV2_,ldaW1_,ldaW2_,ldaM1_,ldaM2_,ldaN1_,ldaN2_,stride_;

  VTest vtest_;
  WTest wtest_;
  MTest mtest_;
  NTest ntest_;

  static void SetUpTestCase()
  {
    TestWithType<_ST_>::SetUpTestCase();
    KernelTestWithMap<_N_>::SetUpTestCase();
    VTest::SetUpTestCase();
    WTest::SetUpTestCase();
    MTest::SetUpTestCase();
    NTest::SetUpTestCase();
  }
  
  /*! Set up routine.
   */
  virtual void SetUp()
  {
    KernelTestWithMap<_N_>::SetUp();
    if (typeImplemented_ && !problemTooSmall_)
    {
      vtest_.SetUp();
      wtest_.SetUp();
      // these probably have a different communicator, but this shouldn't matter here...
      mtest_.SetUp();
      ntest_.SetUp();


      V1_    = vtest_.vec1_;
      V1_vp_ = vtest_.vec1_vp_;
      ldaV1_ = vtest_.lda_;

      V2_    = vtest_.vec2_;
      V2_vp_ = vtest_.vec2_vp_;
      ldaV2_ = vtest_.lda_;


      W1_    = wtest_.vec1_;
      W1_vp_ = wtest_.vec1_vp_;
      ldaW1_ = wtest_.lda_;

      W2_    = wtest_.vec2_;
      W2_vp_ = wtest_.vec2_vp_;
      ldaW2_ = wtest_.lda_;


      M1_    = mtest_.mat1_;
      M1_vp_ = mtest_.mat1_vp_;
      ldaM1_ = mtest_.m_lda_;

      M2_    = mtest_.mat2_;
      M2_vp_ = mtest_.mat2_vp_;
      ldaM2_ = mtest_.m_lda_;


      N1_    = ntest_.mat1_;
      N1_vp_ = ntest_.mat1_vp_;
      ldaN1_ = ntest_.m_lda_;

      N2_    = ntest_.mat2_;
      N2_vp_ = ntest_.mat2_vp_;
      ldaN2_ = ntest_.m_lda_;

    }
    else
    {
      V1_ = NULL;
      V1_vp_ = NULL;
      V2_ = NULL;
      V2_vp_ = NULL;
      W1_ = NULL;
      W1_vp_ = NULL;
      W2_ = NULL;
      W2_vp_ = NULL;
      M1_ = NULL;
      M1_vp_ = NULL;
      M2_ = NULL;
      M2_vp_ = NULL;
      N1_ = NULL;
      N1_vp_ = NULL;
      N2_ = NULL;
      N2_vp_ = NULL;
    }
    stride_=1;
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ntest_.TearDown();
      mtest_.TearDown();
      wtest_.TearDown();
      vtest_.TearDown();
    }
    KernelTestWithMap<_N_>::TearDown();
  }

  static void TearDownTestCase()
  {
    NTest::TearDownTestCase();
    MTest::TearDownTestCase();
    WTest::TearDownTestCase();
    VTest::TearDownTestCase();
    KernelTestWithMap<_N_>::TearDownTestCase();
  }
  
};

  // check augmented kernel with random data
  TEST_F(CLASSNAME, fused_mvsdi_mvTmv) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = ST(2)+st::prand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      SUBR(fused_mvsdi_mvTmv)(alpha,V1_,W1_,N1_,beta,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(alpha,V2_,W2_,beta,M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_NEAR(mt::one(), MvecsEqual(W1_,W2_,mt::one()), sqrt(mt::eps()));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(M1_,M2_), sqrt(mt::eps()));
      
      SUBR(sdMat_from_device)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(M1_,M2_), sqrt(mt::eps()));
    }
  }


#if( _K_ == _M_ )
  // check augmented kernel with random data
  TEST_F(CLASSNAME, fused_mvsdi_mvTmv_self) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V, M, N with random data
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = ST(2)+st::prand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      SUBR(fused_mvsdi_mvTmv)(alpha,V1_,V1_,N1_,beta,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat_inplace)(V2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(alpha,V2_,V2_,beta,M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_NEAR(mt::one(), MvecsEqual(V1_,V2_,mt::one()), 10000*mt::eps());
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(M1_,M2_), sqrt(mt::eps()));

      SUBR(sdMat_from_device)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(M1_,M2_), sqrt(mt::eps()));
    }
  }
#endif


  // check augmented kernel with random data
  TEST_F(CLASSNAME, fused_mvsd_mvTmv) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = ST(2)+st::prand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      SUBR(fused_mvsd_mvTmv)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat)(alpha,V2_,M2_,beta,W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),W2_,W2_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_NEAR(mt::one(), MvecsEqual(W1_,W2_,mt::one()), 1000*mt::eps());
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(N1_,N2_), sqrt(mt::eps()));

        SUBR(sdMat_from_device)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(N1_,N2_), sqrt(mt::eps()));

    }
  }

  // check augmented kernel with result 0
  TEST_F(CLASSNAME, fused_mvsd_mvTmv_result_zero)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      // input data should be ignored by this kernel
      SUBR(sdMat_put_value)(N1_,(_ST_)42.0,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel: W1=alpha*V1*M1+beta*W1 (=0), N1=W1'*W1(=0)
      _ST_ alpha=st::prand(),beta=st::zero();
      SUBR(fused_mvsd_mvTmv)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecEqual(W1_,st::zero()));
      ASSERT_REAL_EQ(mt::one(), SdMatEqual(N1_,st::zero()));
    }
  }

  // check augmented kernel with result 0 but including an actual vector update (beta!=0)
  TEST_F(CLASSNAME, fused_mvsd_mvTmv_update_with_result_zero)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_identity)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,st::zero(),W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // input data should be ignored by this kernel
      SUBR(sdMat_put_value)(N1_,(_ST_)42,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel: W1=alpha*V1(*I)-alpha*W1 (=0), N1=W1'*W1(=0)
      _ST_ alpha=st::prand(),beta=-alpha;
      SUBR(fused_mvsd_mvTmv)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecEqual(W1_,st::zero()));
      ASSERT_REAL_EQ(mt::one(), SdMatEqual(N1_,st::zero()));
    }
  }


  // check augmented kernel with random data
  TEST_F(CLASSNAME, fused_mvsd_mvTmv_nt) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = 0.;

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      SUBR(fused_mvsd_mvTmv)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat)(alpha,V2_,M2_,beta,W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),W2_,W2_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_NEAR(mt::one(), MvecsEqual(W1_,W2_,mt::one()), 10000*mt::eps());
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(M1_,M2_), sqrt(mt::eps()));
    }
  }


  // check augmented kernel with random data
  TEST_F(CLASSNAME, mvec_times_sdMat_add_mvec_times_sdMat) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(V1_,M1_,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),V2_,M2_,st::one(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_NEAR(mt::one(), MvecsEqual(W1_,W2_,mt::one()), mt::sqrt(mt::eps()));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
    }
  }


  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, fused_mvsdi_mvTmv_prec)
#else
  TEST_F(CLASSNAME, DISABLED_fused_mvsdi_mvTmv_prec)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = ST(2)+st::prand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(fused_mvsdi_mvTmv)(alpha,V1_,W1_,N1_,beta,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(alpha,V2_,W2_,beta,M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(W1_,W2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
    }
  }

#if (_M_==_K_)

  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, fused_mvsdi_mvTmv_prec_self)
#else
  TEST_F(CLASSNAME, DISABLED_fused_mvsdi_mvTmv_prec_self)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V, M, N with random data
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = ST(2)+st::prand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(fused_mvsdi_mvTmv)(alpha,V1_,V1_,N1_,beta,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(V2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(alpha,V2_,V2_,beta,M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(V1_vp_,V2_vp_,nloc_,m_,ldaV1_,stride_,vflag_));
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(N1_vp_,N2_vp_,k_,k_,ldaN1_,stride_,mflag_));
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(M1_vp_,M2_vp_,m_,k_,ldaM1_,stride_,mflag_));
    }
  }

#endif

  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, fused_mvsd_mvTmv_prec)
#else
  TEST_F(CLASSNAME, DISABLED_fused_mvsd_mvTmv_prec)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = ST(2)+st::prand(), beta = ST(2)+st::prand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(fused_mvsd_mvTmv)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(alpha,V2_,M2_,beta,W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),W2_,W2_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(W1_,W2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(N1_,N2_), 100*mt::eps());
    }
  }


  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, mvec_times_sdMat_add_mvec_times_sdMat_prec)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_add_mvec_times_sdMat_prec)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(V1_,M1_,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(st::one(),V2_,M2_,st::one(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_NEAR(mt::one(), MvecsEqual(W1_,W2_), sqrt(mt::eps()));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
    }
  }

  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, mvec_times_sdMat_add_mvec_times_sdMat_prec_N0)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_add_mvec_times_sdMat_prec_N0)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(N1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(V1_,M1_,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(st::one(),V2_,M2_,st::one(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(W1_,W2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
    }
  }


  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, mvec_times_sdMat_add_mvec_times_sdMat_prec_M0)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_add_mvec_times_sdMat_prec_M0)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill W, V, M, N with random data
      SUBR(mvec_random)(W1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // copy data
      SUBR(mvec_add_mvec)(st::one(),W1_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(V1_,M1_,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(st::one(),V2_,M2_,st::one(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(W1_,W2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
    }
  }



const int CLASSNAME::k_;
const int CLASSNAME::m_;
const int CLASSNAME::n_;
