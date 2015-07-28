#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithType< _ST_ >,
                 public virtual KernelTestWithMap<_N_>
{

public:

  typedef KernelTestWithVectors<_ST_,_N_,_M_> V1Test;
  typedef KernelTestWithVectors<_ST_,_N_,_K_> V2Test;
  typedef KernelTestWithSdMats<_ST_,_M_,_K_> MTest;

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  static const int k_=_K_;
  
  //! V is n x m
  TYPE(mvec_ptr) V1_,V2_;

  //! W is n x k
  TYPE(mvec_ptr) W1_,W2_;

  //! M is m x k
  TYPE(sdMat_ptr) M1_,M2_;

  //! N is k x k
  TYPE(sdMat_ptr) N1_,N2_;
  
  _ST_ *V1_vp_,*V2_vp_,*W1_vp_,*W2_vp_,*M1_vp_,*M2_vp_,*N1_vp_,*N2_vp_;
  
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  lidx_t ldaV1_,ldaV2_,ldaW1_,ldaW2_,ldaM1_,ldaM2_,ldaN1_,ldaN2_,stride_;
  
  /*! Set up routine.
   */
  virtual void SetUp()
  {
    KernelTestWithType< _ST_ >::SetUp();
    KernelTestWithMap<_N_>::SetUp();
    if (this->typeImplemented_)
    {
      // create vectors V1 and V2, and vector views for setting/checking entries
      PHISTTEST_MVEC_CREATE(&V1_,this->map_,this->m_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_put_value)(V1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_extract_view)(V1_,&V1_vp_,&ldaV1_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&V2_,this->map_,this->m_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_put_value)(V2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_extract_view)(V2_,&V2_vp_,&ldaV2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      // create vectors W1 and W2, and vector views for setting/checking entries
      PHISTTEST_MVEC_CREATE(&W1_,this->map_,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_put_value)(W1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_extract_view)(W1_,&W1_vp_,&ldaW1_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&W2_,this->map_,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_put_value)(W2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_extract_view)(W2_,&W2_vp_,&ldaW2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      // create matrices M1, M2 and views.
      SUBR(sdMat_create)(&M1_,this->m_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(M1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(M1_,&M1_vp_,&this->ldaM1_,&this->iflag_);
      SUBR(sdMat_create)(&M2_,this->m_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(M2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(M2_,&M2_vp_,&this->ldaM2_,&this->iflag_);
      // create matrices N1, N2 and views.
      SUBR(sdMat_create)(&N1_,this->k_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(N1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(N1_,&N1_vp_,&this->ldaN1_,&this->iflag_);
      SUBR(sdMat_create)(&N2_,this->k_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(N2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(N2_,&N2_vp_,&this->ldaN2_,&this->iflag_);

    }
    stride_=1;
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    if (this->typeImplemented_)
    {
      SUBR(mvec_delete)(V1_,&iflag_);
      SUBR(mvec_delete)(V2_,&iflag_);
      SUBR(mvec_delete)(W1_,&iflag_);
      SUBR(mvec_delete)(W2_,&iflag_);
      SUBR(sdMat_delete)(M1_,&iflag_);
      SUBR(sdMat_delete)(M2_,&iflag_);
      SUBR(sdMat_delete)(N1_,&iflag_);
      SUBR(sdMat_delete)(N2_,&iflag_);
    }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
  }

};

  // check augmented kernel with random data
  TEST_F(CLASSNAME, mvecT_times_mvec_times_sdMat_inplace) 
  {
    if (typeImplemented_)
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
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 2.*st::one()+st::rand();

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
      SUBR(mvecT_times_mvec_times_sdMat_inplace)(alpha,V1_,W1_,N1_,beta,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat_inplace)(W2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(alpha,V2_,W2_,beta,M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_NEAR(mt::one(), MvecsEqual(W1_,W2_), sqrt(mt::eps()));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(M1_,M2_), sqrt(mt::eps()));
    }
  }


#if( _K_ == _M_ )
  // check augmented kernel with random data
  TEST_F(CLASSNAME, mvecT_times_mvec_times_sdMat_inplace_self) 
  {
    if (typeImplemented_)
    {
      // fill V, M, N with random data
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 2.*st::one()+st::rand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      SUBR(mvecT_times_mvec_times_sdMat_inplace)(alpha,V1_,V1_,N1_,beta,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat_inplace)(V2_,N2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(alpha,V2_,V2_,beta,M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_NEAR(mt::one(), ArraysEqual(V1_vp_,V2_vp_,nloc_,m_,ldaV1_,stride_,vflag_), mt::eps());
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(N1_vp_,N2_vp_,k_,k_,ldaN1_,stride_,mflag_));
      ASSERT_NEAR(mt::one(), ArraysEqual(M1_vp_,M2_vp_,m_,k_,ldaM1_,stride_,mflag_), sqrt(mt::eps()));
    }
  }
#endif


  // check augmented kernel with random data
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_augmented) 
  {
    if (typeImplemented_)
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
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 2.*st::one()+st::rand();

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
      SUBR(mvec_times_sdMat_augmented)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat)(alpha,V2_,M2_,beta,W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),W2_,W2_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(V1_,V2_));
      ASSERT_REAL_EQ(mt::one(), MvecsEqual(W1_,W2_));
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(M1_,M2_));
      ASSERT_NEAR(mt::one(), SdMatsEqual(N1_,N2_), sqrt(mt::eps()));
    }
  }


  // check augmented kernel with random data
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_augmented_nt) 
  {
    if (typeImplemented_)
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
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 0.;

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
      SUBR(mvec_times_sdMat_augmented)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // do the same calculation by hand
      SUBR(mvec_times_sdMat)(alpha,V2_,M2_,beta,W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),W2_,W2_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // Compare results
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(V1_vp_,V2_vp_,nloc_,m_,ldaV1_,stride_,vflag_));
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(W1_vp_,W2_vp_,nloc_,k_,ldaW1_,stride_,vflag_));
      ASSERT_REAL_EQ(mt::one(), ArraysEqual(M1_vp_,M2_vp_,m_,k_,ldaM1_,stride_,mflag_));
      ASSERT_NEAR(mt::one(), ArraysEqual(N1_vp_,N2_vp_,k_,k_,ldaN1_,stride_,mflag_), sqrt(mt::eps()));
    }
  }


#if( _K_ == 1 || _K_ == 2 || _K_ == 4 )
  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, mvecT_times_mvec_times_sdMat_inplace_prec)
#else
  TEST_F(CLASSNAME, DISABLED_mvecT_times_mvec_times_sdMat_inplace_prec)
#endif
  {
    if (typeImplemented_)
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
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 2.*st::one()+st::rand();

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
      SUBR(mvecT_times_mvec_times_sdMat_inplace)(alpha,V1_,W1_,N1_,beta,M1_,&iflag_);
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


#if( ( _K_ == 1 || _K_ == 2 || _K_ == 4 ) && _K_ == _M_ )
  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, mvecT_times_mvec_times_sdMat_inplace_prec_self)
#else
  TEST_F(CLASSNAME, DISABLED_mvecT_times_mvec_times_sdMat_inplace_prec_self)
#endif
  {
    if (typeImplemented_)
    {
      // fill V, M, N with random data
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(N1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 2.*st::one()+st::rand();

      // copy data
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),M1_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(st::one(),N1_,st::zero(),N2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // actually call augmented kernel
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec_times_sdMat_inplace)(alpha,V1_,V1_,N1_,beta,M1_,&iflag_);
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

#if( _K_ == 1 || _K_ == 2 || _K_ == 4 )
  // check augmented kernel with random data
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, mvec_times_sdMat_augmented_prec)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_augmented_prec)
#endif
  {
    if (typeImplemented_)
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
      _ST_ alpha = 2.*st::one()+st::rand(), beta = 2.*st::one()+st::rand();

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
      SUBR(mvec_times_sdMat_augmented)(alpha,V1_,M1_,beta,W1_,N1_,&iflag_);
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
      ASSERT_REAL_EQ(mt::one(), SdMatsEqual(N1_,N2_));
    }
  }
#endif


#endif /* high-prec test */


const int CLASSNAME::k_;
const int CLASSNAME::m_;
const int CLASSNAME::n_;
