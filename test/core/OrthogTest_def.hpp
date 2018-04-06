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
#error "file not included correctly"
#endif

#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif

#if defined(PHIST_HAVE_BELOS)||defined(PHIST_HAVE_ANASAZI)
#include "phist_trilinos_type_config.h"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual TestWithType< _ST_ >
#ifdef ORTHOG_WITH_HPD_B
                 , public virtual KernelTestWithMassMat<_ST_,_N_>
#endif
                 , public virtual KernelTestWithMap<_N_>
{

public:

  typedef KernelTestWithVectors<_ST_,_N_,_M_> VTest;
  typedef KernelTestWithVectors<_ST_,_N_,_K_> WTest;
#ifdef ORTHOG_WITH_HPD_B
  typedef KernelTestWithMassMat<_ST_,_N_> BTest;
#endif

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  static const int k_=_K_;
  
  //! V is n x m, W is n x k
  TYPE(mvec_ptr) V_, W_,W2_,Q_;
  
  //! linear operator representing a Hermitian positive definite (hpd) matrix B
  //! which defines the inner product in which we want to orthogonalize.
  //! If nullptr it is the identity matrix, this is the case if ORTHOG_WITH_HPD_B is not defined.
  TYPE(linearOp_ptr) B_op;

  //! vectorspaces above, pre-multiplied by B
  TYPE(mvec_ptr) BV_, BW_,BW2_,BQ_;

  //! R0 is m x m, R1 is k x k, R2 is m x k
  TYPE(sdMat_ptr) R0_, R1_, R2_;

  //! for some tests we need to access the raw data of the R matrices and W2 (tmp space)
  _ST_ *R0_vp_,*R1_vp_,*R2_vp_,*W_vp_,*W2_vp_, *BW_vp_, *BW2_vp_;
  phist_lidx ldaR0_, ldaR1_, ldaR2_, ldaW_, ldaW2_, ldaBW_, ldaBW2_, stride_;
  

  static void SetUpTestCase()
  {
    TestWithType<_ST_>::SetUpTestCase();
#ifdef ORTHOG_WITH_HPD_B
    // creates mass matrix B_ without a given context, and
    // re-initializes the KernelTestWithMap by calling SetUpTestCaseWithMap
    BTest::SetUpTestCase(nullptr);
#else
    KernelTestWithMap<_N_>::SetUpTestCase();
#endif
}
  
  //NOTE: we assume stride_=1 here to make the loops simpler,
  //      it seems reasonable to me to do that because mvecs 
  //      are generally stored in column-major order.

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    TestWithType<_ST_>::SetUp();
#ifdef ORTHOG_WITH_HPD_B
    BTest::SetUp();
#else
    KernelTestWithMap<_N_>::SetUp();
#endif
      if(typeImplemented_ && !problemTooSmall_)
      {

        // disable the test because TSQR will not work.
        // This situation occurs if we have a small matrix (_N_=25, say)
        // and many Q vectors (e.g. 10) with multiple MPI procs.
        int globalTooSmall = _N_ < _M_ || _N_ < _K_;
#ifdef PHIST_HAVE_MPI
        int localTooSmall = nloc_ < _M_ || nloc_ < _K_;
        iflag_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        ASSERT_EQ(0,iflag_);
#endif
        problemTooSmall_ = globalTooSmall != 0;
      }


    if (typeImplemented_ && !problemTooSmall_)
    {
      phist_const_map_ptr map=this->map_;

      // create vectors V, W and vector views for setting/checking entries
      PHISTTEST_MVEC_CREATE(&V_,map,this->m_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&W_,map,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&Q_,map,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&W2_,map,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(W_,&W_vp_,&this->ldaW_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(W2_,&W2_vp_,&this->ldaW2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      
      B_op=nullptr; BV_=V_; BW_=W_; BW2_=BW_; BQ_=Q_;

#ifdef ORTHOG_WITH_HPD_B
      B_op=new TYPE(linearOp);
      SUBR(linearOp_wrap_sparseMat)(B_op,B_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // create spaces pre-multiplied by B
      PHISTTEST_MVEC_CREATE(&BV_,map,this->m_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&BW_,map,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&BW2_,map,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&BQ_,map,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

      SUBR(mvec_extract_view)(BW_,&BW_vp_,&this->ldaBW_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(BW2_,&BW2_vp_,&this->ldaBW2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

#endif
      // create matrices R0,R1, R2 and matrix views for setting/checking entries
      SUBR(sdMat_create)(&R0_,this->m_,this->m_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R0_,&R0_vp_,&this->ldaR0_,&this->iflag_);
      SUBR(sdMat_create)(&R1_,this->k_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R1_,&R1_vp_,&this->ldaR1_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_create)(&R2_,this->m_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R2_,&R2_vp_,&this->ldaR2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
    }
    stride_=1;
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      SUBR(mvec_delete)(V_,&iflag_);
      SUBR(mvec_delete)(W_,&iflag_);
      SUBR(mvec_delete)(W2_,&iflag_);
      SUBR(mvec_delete)(Q_,&iflag_);
#ifdef ORTHOG_WITH_HPD_B
    BTest::TearDown();
      SUBR(mvec_delete)(BV_,&iflag_);
      SUBR(mvec_delete)(BW_,&iflag_);
      SUBR(mvec_delete)(BW2_,&iflag_);
      SUBR(mvec_delete)(BQ_,&iflag_);
      delete B_op;
#endif
      SUBR(sdMat_delete)(R0_,&iflag_);
      SUBR(sdMat_delete)(R1_,&iflag_);
      SUBR(sdMat_delete)(R2_,&iflag_);
    }
    KernelTestWithMap<_N_>::TearDown();
    TestWithType<_ST_>::TearDown();
  }

  static void TearDownTestCase()
  {
#ifdef ORTHOG_WITH_HPD_B
    BTest::TearDownTestCase();
#else
    KernelTestWithMap<_N_>::TearDownTestCase();
#endif
    TestWithType< _ST_ >::TearDownTestCase();
  }

  // primary test routine that calls orthog and checks the result vectors
  // for orthogonality etc. Unit tests below just construct various V and W
  // and their rank so that the correct rank detection can be checked. The
  // other arguments are just memory locations that the test routine should
  // use. BV/BW/BQ are not accessed unless ORTHOG_WITH_HPD_B is defined.
  // It is *not* necessary to compute BV and BW beforehand.
  void doOrthogTests(TYPE(mvec_ptr) V, 
                  TYPE(mvec_ptr) W, 
                  TYPE(mvec_ptr) Q,
                  TYPE(mvec_ptr) BV,
                  TYPE(mvec_ptr) BW,
                  TYPE(mvec_ptr) BQ,
                  TYPE(sdMat_ptr) R0,
                  TYPE(sdMat_ptr) R1,
                  TYPE(sdMat_ptr) R2,
                  int expect_iflagV,
                  int expect_iflagVW,
                  int expectedRankV,
                  int expectedRankVW)
  {

      int iflag_in=PHIST_ORTHOG_RANDOMIZE_NULLSPACE;
#ifdef PHIST_HIGH_PRECISION_KERNELS
      // the check using ASSERT_NEAR(mt::one(),...,tol)
      // requires at least eps(1.0) to make sense
      MT tolV=10*mt::eps();
      MT tolW=10*mt::eps();
      iflag_in|=PHIST_ROBUST_REDUCTIONS;
#else
      MT tolV=(MT)10.*VTest::releps(V);
      MT tolW=(MT)10.*WTest::releps(W);
#endif
      // copy Q=W because orthog() works in-place
      SUBR(mvec_add_mvec)(st::one(),W,st::zero(),Q,&iflag_);
      ASSERT_EQ(0,iflag_);

      // orthogonalize the m columns of V. Test that orthog
      // works if the first argument is nullptr.
      int rankVW=-42;
      iflag_=iflag_in;
      SUBR(orthog)(nullptr,V,B_op,R0,nullptr,1,&rankVW,&iflag_);
      if (iflag_!=+2)
      {
        ASSERT_EQ(expect_iflagV,iflag_);
      }
      ASSERT_EQ(expectedRankV,rankVW);
      
      // check for occurences of Inf or NaN
      ASSERT_EQ(0,MvecContainsInfOrNaN(V));
      ASSERT_EQ(0,SdMatContainsInfOrNaN(R0));

#ifdef ORTHOG_WITH_HPD_B
      ASSERT_TRUE(W!=BW);
      // compute BV after orthogonalizing V
      B_op->apply(st::one(),B_op->A,W,st::zero(),BW,&iflag_);
      ASSERT_EQ(0,iflag_);
      B_op->apply(st::one(),B_op->A,V,st::zero(),BV,&iflag_);
      ASSERT_EQ(0,iflag_);

#endif


      // check wether this worked out
      phist_lidx ldaV;
      ST* V_vp=nullptr;
      SUBR(mvec_extract_view)(V,&V_vp,&ldaV,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int nvec_V;
      SUBR(mvec_num_vectors)(V,&nvec_V,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(VTest::nvec_,nvec_V);

      phist_lidx nloc_V;
      SUBR(mvec_my_length)(V,&nloc_V,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_,nloc_V);
      
      SUBR(mvec_from_device)(V,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef ORTHOG_WITH_HPD_B
      // check wether this worked out
      phist_lidx ldaBV;
      ST* BV_vp=nullptr;
      SUBR(mvec_extract_view)(BV,&BV_vp,&ldaBV,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int nvec_BV;
      SUBR(mvec_num_vectors)(BV,&nvec_BV,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(VTest::nvec_,nvec_BV);

      phist_lidx nloc_BV;
      SUBR(mvec_my_length)(BV,&nloc_BV,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_,nloc_BV);
      
      SUBR(mvec_from_device)(BV,&iflag_);
      ASSERT_EQ(0,iflag_);
      EXPECT_NEAR(mt::one(),VTest::ColsAreBNormalized(V_vp,BV_vp,nloc_,ldaV,ldaBV,stride_,mpi_comm_),tolV);
      ASSERT_NEAR(mt::one(),VTest::ColsAreBOrthogonal(V_vp,BV_vp,nloc_,ldaV,ldaBV,stride_,mpi_comm_),tolV);
#else
      EXPECT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp,nloc_,ldaV,stride_,mpi_comm_),tolV);
      ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp,nloc_,ldaV,stride_,mpi_comm_),tolV);
#endif
      int nsteps=2;

      // now orthogonalize W against V. The result should be such that Q*R1=W-V*R2, Q'*Q=I,V'*Q=0
      rankVW=-42;
      iflag_=iflag_in;
      SUBR(orthog)(V,Q,B_op,R1,R2,nsteps,&rankVW,&iflag_);
      ASSERT_EQ(expect_iflagVW,iflag_);

      // check for occurences of Inf or NaN
      ASSERT_EQ(0,MvecContainsInfOrNaN(V));
      ASSERT_EQ(0,MvecContainsInfOrNaN(Q));
      ASSERT_EQ(0,SdMatContainsInfOrNaN(R1));
      ASSERT_EQ(0,SdMatContainsInfOrNaN(R2));

/*
std::cout<<"V=\n";
SUBR(mvec_print)(V,&iflag_);
std::cout<<"W=\n";
SUBR(mvec_print)(W,&iflag_);
std::cout<<"Q=\n";
SUBR(mvec_print)(Q,&iflag_);
std::cout<<"R1=\n";
SUBR(sdMat_print)(R1,&iflag_);
std::cout<<"R2=\n";
SUBR(sdMat_print)(R2,&iflag_);
*/
      ASSERT_EQ(expectedRankVW,rankVW);
      
      // check orthonormality of Q
      phist_lidx ldaQ;
      ST* Q_vp=nullptr;
      SUBR(mvec_extract_view)(Q,&Q_vp,&ldaQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int nvec_Q;
      SUBR(mvec_num_vectors)(Q,&nvec_Q,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(WTest::nvec_,nvec_Q);

      phist_lidx nloc_Q;
      SUBR(mvec_my_length)(V,&nloc_Q,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_,nloc_Q);
      
      SUBR(mvec_from_device)(Q,&iflag_);
      ASSERT_EQ(0,iflag_);

#ifdef ORTHOG_WITH_HPD_B
      // compute BQ after orthogonalizing W against Q
      B_op->apply(st::one(),B_op->A,Q,st::zero(),BQ,&iflag_);
      ASSERT_EQ(0,iflag_);

      phist_lidx ldaBQ;
      ST* BQ_vp=nullptr;
      SUBR(mvec_extract_view)(BQ,&BQ_vp,&ldaBQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int nvec_BQ;
      SUBR(mvec_num_vectors)(BQ,&nvec_BQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(WTest::nvec_,nvec_BQ);

      phist_lidx nloc_BQ;
      SUBR(mvec_my_length)(V,&nloc_BQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_,nloc_BQ);
      
      SUBR(mvec_from_device)(BQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),WTest::ColsAreBNormalized(Q_vp,BQ_vp,nloc_,ldaQ,ldaBQ,stride_,mpi_comm_),tolW);
      ASSERT_NEAR(mt::one(),WTest::ColsAreBOrthogonal(Q_vp,BQ_vp,nloc_,ldaQ,ldaBQ,stride_,mpi_comm_),tolW);
#else
      ASSERT_NEAR(mt::one(),WTest::ColsAreNormalized(Q_vp,nloc_,ldaQ,stride_,mpi_comm_),tolW);
      ASSERT_NEAR(mt::one(),WTest::ColsAreOrthogonal(Q_vp,nloc_,ldaQ,stride_,mpi_comm_),tolW);
#endif

      // check Q and original V are orthogonal to each other, V'Q=0
      TYPE(sdMat_ptr) VtQ=nullptr;
      SUBR(sdMat_create)(&VtQ,nvec_V,nvec_Q,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SdMatOwner<_ST_> _VtQ(VtQ);

      iflag_=iflag_in;
      SUBR(mvecT_times_mvec)(st::one(),V,BQ,st::zero(),VtQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),SdMatEqual(VtQ,st::zero()),tolV);

      // check the decomposition: Q*R1 = W - V*R2 (compute W2=Q*R1+V*R2-W and compare with 0)
      iflag_=iflag_in;
      SUBR(mvec_times_sdMat)(st::one(),Q,R1,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=iflag_in;
      SUBR(mvec_times_sdMat)(st::one(),V,R2,st::one(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=iflag_in;
      SUBR(mvec_add_mvec)(-st::one(),W,st::one(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),ArrayEqual(W2_vp_,nloc_,k_,ldaW2_,stride_,st::zero(),vflag_),tolW);
  }

};

#ifdef ORTHOG_WITH_HPD_B
  TEST_F(CLASSNAME, B_op_is_hpd)
  {
    int iflag_in=0;
#ifdef PHIST_HIGH_PRECISION_KERNELS
    iflag_in=PHIST_ROBUST_REDUCTIONS;
#endif
    if (problemTooSmall_ || !typeImplemented_) return;
    SUBR(mvec_random)(V_,&iflag_);
    ASSERT_EQ(0,iflag_);
    B_op->apply(st::one(),B_op->A,V_,st::zero(),BV_,&iflag_);
    ASSERT_EQ(0,iflag_);
    iflag_=iflag_in;
    SUBR(mvecT_times_mvec)(st::one(),V_,BV_,st::zero(),R0_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _MT_ sym_err=mt::zero();
    _MT_ min_diag=(_MT_)1.0e10;
    for (int i=0; i<m_; i++)
    {
      min_diag=std::min(min_diag,st::real(R0_vp_[i*ldaR0_+i]));
      for (int j=0; j<m_; j++)
      {
        sym_err = std::max(sym_err, std::abs(st::conj(R0_vp_[i*ldaR0_+j])-R0_vp_[j*ldaR0_+i]));
      }
    }
#ifdef PHIST_HIGH_PRECISION_KERNELS
    ASSERT_REAL_EQ(mt::one(), mt::one()+sym_err);
#else
    ASSERT_NEAR(mt::one(), mt::one()+sym_err,100.0*mt::eps());
#endif
    ASSERT_GT(min_diag,std::sqrt(mt::eps()));
  }
  
  TEST_F(CLASSNAME,orthog_impl_updates_BW)
  {
    if (!problemTooSmall_ && typeImplemented_)
    {
      // create some V and BV, to make it a bit more interesting, make the first and last column linearly dependent.
      SUBR(mvec_random)(V_,&iflag_);
      ASSERT_EQ(0,iflag_);
      int m=VTest::nvec_;
      if (m>1)
      {
        TYPE(mvec_ptr) v1=nullptr, vm=nullptr;
        SUBR(mvec_view_block)(V_,&v1,0,0,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(V_,&vm,m-1,m-1,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(mvec_add_mvec)((_ST_)2,v1,st::zero(),vm,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(mvec_delete)(v1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(vm,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      SUBR(linearOp_apply)(st::one(),B_op,V_,st::zero(),BV_,&iflag_);
      int numSweeps=2;
      _MT_ orthoEps=std::sqrt(mt::eps());
      _MT_ rankTol=mt::eps();
      int rankV;
      TYPE(sdMat_ptr) VtV = nullptr;
      SUBR(sdMat_create)(&VtV,m,m,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SdMatOwner< _ST_ > _VtV(VtV);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V_,BV_,st::zero(),VtV,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(orthog_impl)(nullptr,V_,B_op,BV_,VtV,R0_,nullptr,numSweeps,&rankV,rankTol,orthoEps,&iflag_);
      ASSERT_EQ(m>1?1:0,iflag_); // +1 means not full rank
      ASSERT_EQ(m>1?m-1:1,rankV); // rank is k-1 if there are at least two cols, see above.
      // now check that BV=B*V still holds
      SUBR(linearOp_apply)(-st::one(),B_op,V_,st::one(),BV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),VTest::MvecEqual(BV_,st::zero()),VTest::releps());
    }
  }
#endif

  // check if vectors are normalized correctly after QR factorization
  TEST_F(CLASSNAME, test_with_random_vectors)
  {
    if (typeImplemented_ && !problemTooSmall_)
      {
      // fill V and W with random numbers
      SUBR(mvec_random)(V_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(W_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // test orthog routine, expect full rank of V and [V W]
      doOrthogTests(V_, W_, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        0, 0, m_, m_+k_);
      
    }
  }
  // check if random orthogonal vectors are generated automatically if filled with one-vectors
  TEST_F(CLASSNAME, test_with_one_vectors)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      // fill V and W with one-vectors
      SUBR(mvec_put_value)(V_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(W_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);

      // test orthog routine, expect V and [V, W] to have rank 1
      doOrthogTests(V_, W_, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        m_>1? 1:0, 1, 1, m_);
           
    }
  }

  // check if random orthogonal vectors are generated automatically if filled with one-vectors
  TEST_F(CLASSNAME, test_with_constant_vectors)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      _ST_ v1= ST(0.9344646) - ST(1.04357)*st::cmplx_I();
      _ST_ v2=-ST(0.9554373) + ST(1.08158)*st::cmplx_I();
      // fill V and W with constant entries!=1
      SUBR(mvec_put_value)(V_,v1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(W_,v2,&iflag_);
      ASSERT_EQ(0,iflag_);

      // test orthog routine, expect V and [V, W] to have rank 1
      doOrthogTests(V_, W_, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        m_>1? 1:0, 1, 1, m_);
           
    }
  }


  // check if we can orthogonalize "into a view", that is for instance
  // compute R1, R2 in [R0 R2;
  //                    0  R1]
  TEST_F(CLASSNAME, random_vectors_into_viewed_R)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with random numbers
      SUBR(mvec_random)(V_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(W_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      phist_lidx nR=_M_+_K_;
      TYPE(sdMat_ptr) R=nullptr;
      SUBR(sdMat_create)(&R,nR,nR,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(R0_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R0_=nullptr;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=nullptr;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=nullptr;
      SUBR(sdMat_view_block)(R,&R0_,0,_M_-1,0,_M_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      // R1 is the k x k diagonal block:
      SUBR(sdMat_view_block)(R,&R1_,_M_,_M_+_K_-1,_M_,_M_+_K_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      // R2 is the m'th block col with k cols and m rows
      SUBR(sdMat_view_block)(R,&R2_,0,_M_-1,_M_,_M_+_K_-1,&iflag_);
      ASSERT_EQ(0,iflag_);

      // renew raw views
      SUBR(sdMat_extract_view)(R0_,&R0_vp_,&this->ldaR0_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R1_,&R1_vp_,&this->ldaR1_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R2_,&R2_vp_,&this->ldaR2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

      // test orthog routine, expect full rank of V and [V W]
      doOrthogTests(V_, W_, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        0, 0, m_, m_+k_);

      // check the location of the resulting R parts
      _ST_ *R_vp;
      phist_lidx ldaR;

      SUBR(sdMat_extract_view)(R,&R_vp,&ldaR,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

#ifdef PHIST_SDMATS_ROW_MAJOR
#warning "test not implemented for row-major sdMats"
return;
#endif

        _MT_ errR2=mt::one();
        for (phist_lidx i=0; i<m_; i++)
        {
          for (phist_lidx j=0; j<k_; j++)
          {
            phist_lidx jj = m_+j;
            _ST_ a = R_vp[jj*ldaR+i];
            _ST_ b = R2_vp_[j*ldaR2_+i];
            _MT_ denom=st::abs(a+b);
            if (denom==mt::zero()) denom=mt::one();
            errR2+= st::abs(a-b)/denom;
          }
        }
        _MT_ errR1=mt::one();
        for (phist_lidx i=0; i<k_; i++)
        {
          for (phist_lidx j=0; j<k_; j++)
          {
            phist_lidx ii = m_+i;
            phist_lidx jj = m_+j;
            _ST_ a = R_vp[jj*ldaR+ii];
            _ST_ b = R1_vp_[j*ldaR1_+i];
            _MT_ denom=st::abs(a+b);
            if (denom==mt::zero()) denom=mt::one();
            errR1+= st::abs(a-b)/denom;
          }
        }
#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_DEB("matrix H (last block-col should contain orthog-coefficients)");
      SUBR(sdMat_print)(R,&iflag_);
#endif
      ASSERT_REAL_EQ(mt::one(),errR1);
      ASSERT_REAL_EQ(mt::one(),errR2);

      SUBR(sdMat_delete)(R0_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R0_=nullptr;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=nullptr;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=nullptr;
      SUBR(sdMat_delete)(R,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check if we can orthogonalize "into a view", that is for instance
  // compute R1, R2 in [R0 R2;
  //                    0  R1]
  TEST_F(CLASSNAME, one_vectors_into_viewed_R)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with random numbers
      SUBR(mvec_put_value)(V_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(W_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      
      phist_lidx nR=_M_+_K_;
      TYPE(sdMat_ptr) R=nullptr;
      SUBR(sdMat_create)(&R,nR,nR,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(R0_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R0_=nullptr;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=nullptr;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=nullptr;
      SUBR(sdMat_view_block)(R,&R0_,0,_M_-1,0,_M_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      // R1 is the k x k diagonal block:
      SUBR(sdMat_view_block)(R,&R1_,_M_,_M_+_K_-1,_M_,_M_+_K_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      // R2 is the m'th block col with k cols and m rows
      SUBR(sdMat_view_block)(R,&R2_,0,_M_-1,_M_,_M_+_K_-1,&iflag_);
      ASSERT_EQ(0,iflag_);

      // renew raw views
      SUBR(sdMat_extract_view)(R0_,&R0_vp_,&this->ldaR0_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R1_,&R1_vp_,&this->ldaR1_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_extract_view)(R2_,&R2_vp_,&this->ldaR2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

      // test orthog routine, expect V and [V, W] to have rank 1
      doOrthogTests(V_, W_, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        m_>1? 1:0, 1, 1, m_);
      
      // check the location of the resulting R parts
      _ST_ *R_vp;
      phist_lidx ldaR;

      SUBR(sdMat_extract_view)(R,&R_vp,&ldaR,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

#ifdef PHIST_SDMATS_ROW_MAJOR
#warning "test not implemented for row-major sdMats"
        return;
#endif

        _MT_ errR2=mt::one();
        for (phist_lidx i=0; i<m_; i++)
        {
          for (phist_lidx j=0; j<k_; j++)
          {
            phist_lidx jj = m_+j;
            _ST_ a = R_vp[jj*ldaR+i];
            _ST_ b = R2_vp_[j*ldaR2_+i];
            _MT_ denom=st::abs(a+b);
            if (denom==mt::zero()) denom=mt::one();
            errR2+= st::abs(a-b)/denom;
          }
        }
        _MT_ errR1=mt::one();
        for (phist_lidx i=0; i<k_; i++)
        {
          for (phist_lidx j=0; j<k_; j++)
          {
            phist_lidx ii = m_+i;
            phist_lidx jj = m_+j;
            _ST_ a = R_vp[jj*ldaR+ii];
            _ST_ b = R1_vp_[j*ldaR1_+i];
            _MT_ denom=st::abs(a+b);
            if (denom==mt::zero()) denom=mt::one();
            errR1+= st::abs(a-b)/denom;
          }
        }
#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_DEB("matrix H (last block-col should contain orthog-coefficients)");
      SUBR(sdMat_print)(R,&iflag_);
#endif    
      ASSERT_REAL_EQ(mt::one(),errR1);
      ASSERT_REAL_EQ(mt::one(),errR2);

      SUBR(sdMat_delete)(R0_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R0_=nullptr;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=nullptr;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=nullptr;
      SUBR(sdMat_delete)(R,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check if we can orthogonalize W against V where V and W are views into Z=[V W].
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, DISABLED_random_vectors_viewing_same_block)
#else
  TEST_F(CLASSNAME, random_vectors_viewing_same_block)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      phist_const_map_ptr map=map_;
#ifdef ORTHOG_WITH_HPD_B
      map=B_op->domain_map;
#endif      
      TYPE(mvec_ptr) V_big=nullptr, V=nullptr, W=nullptr;
      int ncols = m_+k_+13;// some extra padding to make the test more interesting
      PHISTTEST_MVEC_CREATE(&V_big,map,ncols,&iflag_);
      ASSERT_EQ(0,iflag_);

      // fill entire block with random numbers
      SUBR(mvec_random)(V_big,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // create views
      SUBR(mvec_view_block)(V_big,&V,1,m_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_view_block)(V_big,&W,m_+1,m_+k_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // test orthog routine, expect full rank of V and [V W]
      doOrthogTests(V, W, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        0, 0, m_, m_+k_);
      
      SUBR(mvec_delete)(V,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(W,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(V_big,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check if we can orthogonalize W against V where V and W are views into Z=[V W].
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, DISABLED_one_vectors_viewing_same_block)
#else
  TEST_F(CLASSNAME, one_vectors_viewing_same_block)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      phist_const_map_ptr map=map_;
#ifdef ORTHOG_WITH_HPD_B
      map=B_op->domain_map;
#endif      
      TYPE(mvec_ptr) V_big=nullptr, V=nullptr, W=nullptr;
      int ncols = m_+k_+13;// some extra padding to make the test more interesting
      PHISTTEST_MVEC_CREATE(&V_big,map,ncols,&iflag_);
      ASSERT_EQ(0,iflag_);

      // fill entire block with ones
      SUBR(mvec_put_value)(V_big,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // create views
      SUBR(mvec_view_block)(V_big,&V,1,m_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_view_block)(V_big,&W,m_+1,m_+k_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // test orthog routine, expect V and [V, W] to have rank 1
      doOrthogTests(V, W, Q_, BV_, BW_, BQ_, R0_, R1_, R2_, 
        m_>1? 1:0, 1, 1, m_);

      SUBR(mvec_delete)(V,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(W,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(V_big,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

#if defined(PHIST_HAVE_BELOS)&&defined(PHIST_TRILINOS_TYPE_AVAIL)

  // compare our overloaded ICGS (based on orthog) with another Belos OrthoManager class.
  TEST_F(CLASSNAME,belos_ortho_manager)
  {
    typedef phist::BelosMV< _ST_ > MV;
    typedef Teuchos::SerialDenseMatrix<int, _ST_ > SDM;
    typedef phist::ScalarTraits< _ST_ >::linearOp_t OP;
    typedef Belos::MultiVecTraits<_ST_, MV> MVT;
    typedef Belos::OperatorTraits<_ST_, MV, OP> OPT;
  
    if (typeImplemented_ && !problemTooSmall_)
    {
      // create a random orthogonal space V with m columns
      SUBR(mvec_random)(V_,&iflag_);
      iflag_=PHIST_ORTHOG_RANDOMIZE_NULLSPACE;
      int rankV0;
      SUBR(orthog)(nullptr,V_,B_op,R0_,nullptr,1,&rankV0,&iflag_);
      ASSERT_EQ(0,iflag_);
 
      // random vector block W with k columns, backup to W2.
      SUBR(mvec_random)(W_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),W_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // copy Q_=W_ and use it as a backup. Note the different notations:
      // We use V as the orthogonal basis, in Belos it's called Q, and our W is called V in Belos.
      // What we call R1 and R2 are called B and C, respectively.
      // So after projectAndNormalize, W_orig-V*C=W*B
      SUBR(mvec_add_mvec)(st::one(),W_,st::zero(),Q_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),BW_,st::zero(),BW2_,&iflag_);
      ASSERT_EQ(0,iflag_);

#ifdef ORTHOG_WITH_HPD_B
      // compute BV after orthogonalizing V
      B_op->apply(st::one(),B_op->A,W_,st::zero(),BW_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      // create the Belos objects
      Teuchos::RCP<const MV> V =phist::mvec_rcp< _ST_ >((TYPE(const_mvec_ptr))V_, false);
      Teuchos::RCP<MV> W=phist::mvec_rcp< _ST_ >( W_,false),
                       W0=phist::mvec_rcp< _ST_ >(Q_,false),
                       W2=phist::mvec_rcp< _ST_ >(W2_,false), 
                       BW=phist::mvec_rcp< _ST_ >(BW_,false),
                       BW2=phist::mvec_rcp< _ST_ >(BW2_,false);
      Teuchos::RCP<OP> Op=Teuchos::rcp(B_op,false);
     
      Teuchos::RCP<SDM> B = Teuchos::rcp(new SDM(k_,k_)), B2=Teuchos::rcp(new SDM(k_,k_));
      Teuchos::RCP<SDM> C = Teuchos::rcp(new SDM(m_,k_)), C2=Teuchos::rcp(new SDM(m_,k_));
      
      Teuchos::Array<Teuchos::RCP<SDM> > C_array(1,C),C2_array(1,C2);
      Teuchos::ArrayView<Teuchos::RCP<const MV> > V_array(&V,1);
      
      int max_blk_ortho=2;
      _MT_ blk_tol=std::sqrt(mt::eps());
      _MT_ dep_tol=mt::eps();
      _MT_ sing_tol=mt::eps();
      
      Belos::ICGSOrthoManager<_ST_,MV,OP> myOrtho("phist/orthog",Op,
                                                  max_blk_ortho, blk_tol, sing_tol);
/* note: it seems that this ortho-manager doesn't treat the operator Op correctly
      Belos::DGKSOrthoManager<_ST_,MV,OP> theirOrtho("Belos/DGKS",Op,
                                                  max_blk_ortho, blk_tol, dep_tol, sing_tol);
*/
      ::Belos::IMGSOrthoManager<_ST_,MV,OP> theirOrtho("Belos/IMGS",Op,
                                                  max_blk_ortho, blk_tol, sing_tol);

      int my_ret=myOrtho.projectAndNormalize(*W,BW,C_array,B,V_array);
      int their_ret=theirOrtho.projectAndNormalize(*W2,BW2,C2_array,B2,V_array);

      // should return the same code
      ASSERT_EQ(my_ret,their_ret);

      std::cout << "my B = "<<*B << std::endl;
      std::cout << "their B = "<<*B2 << std::endl;
      std::cout << "my C = "<<*C << std::endl;
      std::cout << "their C = "<<*C2 << std::endl;

      // there may be sign switches in B, depending on the algorithm used, so 'normalize' the signs first.
      std::vector<_ST_> col_signs(k_),col_signs2(k_);
      for (int i=0; i<k_; i++)
      {
        col_signs[i] = (_ST_)((mt::zero() < st::real((*B)(i,i))) - (st::real((*B)(i,i)) < mt::zero()));
        col_signs2[i] = (_ST_)((mt::zero() < st::real((*B2)(i,i))) - (st::real((*B2)(i,i)) < mt::zero()));
        for (int j=0; j<B->numCols(); j++)
        {
          (*B)(i,j)*=col_signs[i];
          (*B2)(i,j)*=col_signs2[i];
        }
      }
      MVT::MvScale(*W,col_signs);
      MVT::MvScale(*W2,col_signs2);
#ifdef ORTHOG_WITH_HPD_B
      MVT::MvScale(*BW,col_signs);
      MVT::MvScale(*BW2,col_signs2);
#endif
      /* compute the orthog relation QB - (W-VC), should be almost zero for both methods */
      Teuchos::RCP<MV> orthogCond = MVT::CloneCopy(*W0);
      Teuchos::RCP<MV> orthogCond2 = MVT::CloneCopy(*W0);
      
      MVT::MvTimesMatAddMv(st::one(), *W, *B,  -st::one(), *orthogCond);
      MVT::MvTimesMatAddMv(st::one(), *W2,*B2, -st::one(), *orthogCond2);

      MVT::MvTimesMatAddMv(st::one(),*V,*C,  st::one(), *orthogCond);
      MVT::MvTimesMatAddMv(st::one(),*V,*C2, st::one(), *orthogCond2);
      
      std::vector< _MT_ > normsOrthogCond(k_),normsOrthogCond2(k_);
      MVT::MvNorm(*orthogCond,normsOrthogCond);
      MVT::MvNorm(*orthogCond2,normsOrthogCond2);
      
      auto maxOrthogCondError=std::max_element(normsOrthogCond.begin(),  normsOrthogCond.end());
      auto maxOrthogCondError2=std::max_element(normsOrthogCond2.begin(),  normsOrthogCond2.end());

      ASSERT_NEAR(mt::one(),mt::one()+*maxOrthogCondError, std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),mt::one()+*maxOrthogCondError2, std::sqrt(mt::eps()));

#ifdef ORTHOG_WITH_HPD_B
      // check that BW after the orthogonalization is correctly updated
      Teuchos::RCP<MV> BW_expl = MVT::Clone(*W,k_);
      Teuchos::RCP<MV> BW_expl2 = MVT::Clone(*W,k_);
      Teuchos::RCP<MV> BW_updated = MVT::Clone(*W,k_);
      Teuchos::RCP<MV> BW_updated2 = MVT::Clone(*W,k_);
      
      OPT::Apply(*Op,*W,*BW_expl,Belos::NOTRANS);
      OPT::Apply(*Op,*W2,*BW_expl2,Belos::NOTRANS);
      
      MVT::MvAddMv(-st::one(),*BW, st::one(),*BW_expl, *BW_updated);
      MVT::MvAddMv(-st::one(),*BW2,st::one(),*BW_expl2, *BW_updated2);

      std::vector< _MT_ > normsBWupdated(k_),normsBWupdated2(k_);
      MVT::MvNorm(*BW_updated,  normsBWupdated);
      MVT::MvNorm(*BW_updated2, normsBWupdated2);
      
      auto maxBWupdatedError=std::max_element(normsBWupdated.begin(),  normsBWupdated.end());
      auto maxBWupdatedError2=std::max_element(normsBWupdated2.begin(),  normsBWupdated2.end());

      ASSERT_NEAR(mt::one(),mt::one()+*maxBWupdatedError, std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),mt::one()+*maxBWupdatedError2, std::sqrt(mt::eps()));
      
#endif      
      // check that W is correctly orthonormalized 
#ifdef ORTHOG_WITH_HPD_B
      ASSERT_NEAR(mt::one(),WTest::ColsAreBOrthogonal(W_vp_,BW_vp_,nloc_,ldaW_,ldaBW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreBNormalized(W_vp_,BW_vp_,nloc_,ldaW_,ldaBW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreBOrthogonal(W2_vp_,BW2_vp_,nloc_,ldaW2_,ldaBW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreBNormalized(W2_vp_,BW2_vp_,nloc_,ldaW2_,ldaBW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
#else
      ASSERT_NEAR(mt::one(),WTest::ColsAreOrthogonal(W_vp_,nloc_,ldaW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreNormalized(W_vp_,nloc_,ldaW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreOrthogonal(W2_vp_,nloc_,ldaW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreNormalized(W2_vp_,nloc_,ldaW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
#endif
      
      SDM Bdiff=*B; Bdiff-=*B2;
      SDM Cdiff=*C; Cdiff-=*C2;

      //std::cout << "Bdiff = "<<Bdiff << std::endl;
      //std::cout << "Cdiff = "<<Cdiff << std::endl;
       
      ASSERT_NEAR(mt::one(), mt::one()+Bdiff.normFrobenius(),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(), mt::one()+Cdiff.normFrobenius(),std::sqrt(mt::eps()));
    }
  }

#endif

// todo!
#if defined(PHIST_HAVE_ANASAZI)&&defined(PHIST_TRILINOS_TYPE_AVAIL)

  // compare our overloaded SVQB (based on orthog) with another Anasazi OrthoManager class.
  TEST_F(CLASSNAME,anasazi_ortho_manager)
  {
    typedef phist::BelosMV< _ST_ > MV;
    typedef Teuchos::SerialDenseMatrix<int, _ST_ > SDM;
    typedef phist::ScalarTraits< _ST_ >::linearOp_t OP;
    typedef Belos::MultiVecTraits<_ST_, MV> MVT;
    typedef Belos::OperatorTraits<_ST_, MV, OP> OPT;
  
    if (typeImplemented_ && !problemTooSmall_)
    {
      // create a random orthogonal space V with m columns
      SUBR(mvec_random)(V_,&iflag_);
      iflag_=PHIST_ORTHOG_RANDOMIZE_NULLSPACE;
      int rankV0;
      SUBR(orthog)(nullptr,V_,B_op,R0_,nullptr,1,&rankV0,&iflag_);
      ASSERT_EQ(0,iflag_);
 
      // random vector block W with k columns, backup to W2.
      SUBR(mvec_random)(W_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),W_,st::zero(),W2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // copy Q_=W_ and use it as a backup. Note the different notations:
      // We use V as the orthogonal basis, in Belos it's called Q, and our W is called V in Belos.
      // What we call R1 and R2 are called B and C, respectively.
      // So after projectAndNormalize, W_orig-V*C=W*B
      SUBR(mvec_add_mvec)(st::one(),W_,st::zero(),Q_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),BW_,st::zero(),BW2_,&iflag_);
      ASSERT_EQ(0,iflag_);

#ifdef ORTHOG_WITH_HPD_B
      // compute BV after orthogonalizing V
      B_op->apply(st::one(),B_op->A,W_,st::zero(),BW_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      // create the Belos objects
      Teuchos::RCP<const MV> V =phist::mvec_rcp< _ST_ >((TYPE(const_mvec_ptr))V_, false);
      Teuchos::RCP<MV> W=phist::mvec_rcp< _ST_ >( W_,false),
                       W0=phist::mvec_rcp< _ST_ >(Q_,false),
                       W2=phist::mvec_rcp< _ST_ >(W2_,false), 
                       BW=phist::mvec_rcp< _ST_ >(BW_,false),
                       BW2=phist::mvec_rcp< _ST_ >(BW2_,false);
      Teuchos::RCP<OP> Op=Teuchos::rcp(B_op,false);
     
      Teuchos::RCP<SDM> B = Teuchos::rcp(new SDM(k_,k_)), B2=Teuchos::rcp(new SDM(k_,k_));
      Teuchos::RCP<SDM> C = Teuchos::rcp(new SDM(m_,k_)), C2=Teuchos::rcp(new SDM(m_,k_));
      
      Teuchos::Array<Teuchos::RCP<SDM> > C_array(1,C),C2_array(1,C2);
      Teuchos::Array<Teuchos::RCP<const MV> > V_array(1,V);
      
      int max_blk_ortho=2;
      _MT_ blk_tol=std::sqrt(mt::eps());
      _MT_ dep_tol=mt::eps();
      _MT_ sing_tol=mt::eps();

      Anasazi::SVQBOrthoManager<_ST_,MV,OP> myOrtho(Op,false);
      ::Anasazi::ICGSOrthoManager<_ST_,MV,OP> theirOrtho(Op,max_blk_ortho, sing_tol, blk_tol);

      int my_ret=myOrtho.projectAndNormalize(*W,V_array,C_array,B);
      int their_ret=theirOrtho.projectAndNormalize(*W2,V_array,C2_array,B2);
      
      // should return the same code
      ASSERT_EQ(my_ret,their_ret);

      std::cout << "my B = "<<*B << std::endl;
      std::cout << "their B = "<<*B2 << std::endl;
      std::cout << "my C = "<<*C << std::endl;
      std::cout << "their C = "<<*C2 << std::endl;

      // there may be sign switches in B, depending on the algorithm used, so 'normalize' the signs first.
      std::vector<_ST_> col_signs(k_),col_signs2(k_);
      for (int i=0; i<k_; i++)
      {
        col_signs[i] = (_ST_)((mt::zero() < st::real((*B)(i,i))) - (st::real((*B)(i,i)) < mt::zero()));
        col_signs2[i] = (_ST_)((mt::zero() < st::real((*B2)(i,i))) - (st::real((*B2)(i,i)) < mt::zero()));
        for (int j=0; j<B->numCols(); j++)
        {
          (*B)(i,j)*=col_signs[i];
          (*B2)(i,j)*=col_signs2[i];
        }
      }
      MVT::MvScale(*W,col_signs);
      MVT::MvScale(*W2,col_signs2);
#ifdef ORTHOG_WITH_HPD_B
      MVT::MvScale(*BW,col_signs);
      MVT::MvScale(*BW2,col_signs2);
#endif
      /* compute the orthog relation QB - (W-VC), should be almost zero for both methods */
      Teuchos::RCP<MV> orthogCond = MVT::CloneCopy(*W0);
      Teuchos::RCP<MV> orthogCond2 = MVT::CloneCopy(*W0);
      
      MVT::MvTimesMatAddMv(st::one(), *W, *B,  -st::one(), *orthogCond);
      MVT::MvTimesMatAddMv(st::one(), *W2,*B2, -st::one(), *orthogCond2);

      MVT::MvTimesMatAddMv(st::one(),*V,*C,  st::one(), *orthogCond);
      MVT::MvTimesMatAddMv(st::one(),*V,*C2, st::one(), *orthogCond2);
      
      std::vector< _MT_ > normsOrthogCond(k_),normsOrthogCond2(k_);
      MVT::MvNorm(*orthogCond,normsOrthogCond);
      MVT::MvNorm(*orthogCond2,normsOrthogCond2);
      
      auto maxOrthogCondError=std::max_element(normsOrthogCond.begin(),  normsOrthogCond.end());
      auto maxOrthogCondError2=std::max_element(normsOrthogCond2.begin(),  normsOrthogCond2.end());

      ASSERT_NEAR(mt::one(),mt::one()+*maxOrthogCondError, std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),mt::one()+*maxOrthogCondError2, std::sqrt(mt::eps()));

#ifdef ORTHOG_WITH_HPD_B
      // check that BW after the orthogonalization is correctly updated
      Teuchos::RCP<MV> BW_expl = MVT::Clone(*W,k_);
      Teuchos::RCP<MV> BW_expl2 = MVT::Clone(*W,k_);
      Teuchos::RCP<MV> BW_updated = MVT::Clone(*W,k_);
      Teuchos::RCP<MV> BW_updated2 = MVT::Clone(*W,k_);
      
      OPT::Apply(*Op,*W,*BW_expl,Belos::NOTRANS);
      OPT::Apply(*Op,*W2,*BW_expl2,Belos::NOTRANS);
      
      MVT::MvAddMv(-st::one(),*BW, st::one(),*BW_expl, *BW_updated);
      MVT::MvAddMv(-st::one(),*BW2,st::one(),*BW_expl2, *BW_updated2);

      std::vector< _MT_ > normsBWupdated(k_),normsBWupdated2(k_);
      MVT::MvNorm(*BW_updated,  normsBWupdated);
      MVT::MvNorm(*BW_updated2, normsBWupdated2);
      
      auto maxBWupdatedError=std::max_element(normsBWupdated.begin(),  normsBWupdated.end());
      auto maxBWupdatedError2=std::max_element(normsBWupdated2.begin(),  normsBWupdated2.end());

      ASSERT_NEAR(mt::one(),mt::one()+*maxBWupdatedError, std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),mt::one()+*maxBWupdatedError2, std::sqrt(mt::eps()));
      
#endif      
      // check that W is correctly orthonormalized 
#ifdef ORTHOG_WITH_HPD_B
      ASSERT_NEAR(mt::one(),WTest::ColsAreBOrthogonal(W_vp_,BW_vp_,nloc_,ldaW_,ldaBW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreBNormalized(W_vp_,BW_vp_,nloc_,ldaW_,ldaBW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreBOrthogonal(W2_vp_,BW2_vp_,nloc_,ldaW2_,ldaBW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreBNormalized(W2_vp_,BW2_vp_,nloc_,ldaW2_,ldaBW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
#else
      ASSERT_NEAR(mt::one(),WTest::ColsAreOrthogonal(W_vp_,nloc_,ldaW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreNormalized(W_vp_,nloc_,ldaW_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreOrthogonal(W2_vp_,nloc_,ldaW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(),WTest::ColsAreNormalized(W2_vp_,nloc_,ldaW2_,stride_,mpi_comm_),std::sqrt(mt::eps()));
#endif
      
      SDM Bdiff=*B; Bdiff-=*B2;
      SDM Cdiff=*C; Cdiff-=*C2;

      //std::cout << "Bdiff = "<<Bdiff << std::endl;
      //std::cout << "Cdiff = "<<Cdiff << std::endl;
       
      ASSERT_NEAR(mt::one(), mt::one()+Bdiff.normFrobenius(),std::sqrt(mt::eps()));
      ASSERT_NEAR(mt::one(), mt::one()+Cdiff.normFrobenius(),std::sqrt(mt::eps()));
    }
  }

#endif
