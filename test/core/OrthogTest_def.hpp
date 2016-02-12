#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual TestWithType< _ST_ >,
                 public virtual KernelTestWithMap<_N_>
  {

public:

  typedef KernelTestWithVectors<_ST_,_N_,_M_> VTest;
  typedef KernelTestWithVectors<_ST_,_N_,_K_> WTest;

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  static const int k_=_K_;
  
  //! V is n x m, W is n x k  
  TYPE(mvec_ptr) V_, W_,W2_,Q_;

  //! R1 is m x m, R2 is k x k
  TYPE(sdMat_ptr) R0_, R1_, R2_;
  
  _ST_ *V_vp_,*W_vp_,*W2_vp_,*Q_vp_,*R0_vp_,*R1_vp_,*R2_vp_;
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  lidx_t ldaV_,ldaW_,ldaW2_,ldaQ_,ldaR0_,ldaR1_,ldaR2_,stride_;


  static void SetUpTestCase()
  {
    KernelTestWithMap<_N_>::SetUpTestCase();
    TestWithType<_ST_>::SetUpTestCase();
  }
  
  //NOTE: we assume stride_=1 here to make the loops simpler,
  //      it seems reasonable to me to do that because mvecs 
  //      are generally stored in column-major order.

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    TestWithType<_ST_>::SetUp();
    KernelTestWithMap<_N_>::SetUp();

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
      // create vectors V, W and vector views for setting/checking entries
      PHISTTEST_MVEC_CREATE(&V_,this->map_,this->m_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&W_,this->map_,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(W_,&W_vp_,&ldaW_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&Q_,this->map_,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(Q_,&Q_vp_,&ldaQ_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      PHISTTEST_MVEC_CREATE(&W2_,this->map_,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_extract_view)(W2_,&W2_vp_,&ldaW2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
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
      SUBR(sdMat_delete)(R0_,&iflag_);
      SUBR(sdMat_delete)(R1_,&iflag_);
      SUBR(sdMat_delete)(R2_,&iflag_);
      }
    KernelTestWithMap<_N_>::TearDown();
    TestWithType<_ST_>::TearDown();
    }

  void doOrthogTests(TYPE(mvec_ptr) V, 
                  TYPE(mvec_ptr) W, 
                  TYPE(mvec_ptr) Q,
                  TYPE(sdMat_ptr) R0,
                  TYPE(sdMat_ptr) R1,
                  TYPE(sdMat_ptr) R2,
                  int expect_iflagV,
                  int expect_iflagVW,
                  int expectedRankV,
                  int expectedRankVW)
  {

      MT tolV=(MT)10.*VTest::releps(V);
      MT tolW=(MT)10.*WTest::releps(W);

      // copy Q=W because orthog() works in-place
      SUBR(mvec_add_mvec)(st::one(),W,st::zero(),Q,&iflag_);
      ASSERT_EQ(0,iflag_);
      // orthogonalize the m columns of V. Test that orthog
      // works if the first argument is NULL.
      int rankVW=-42;
      SUBR(orthog)(NULL,V,NULL,R0,NULL,1,&rankVW,&iflag_);
      if (iflag_!=+2)
      {
        ASSERT_EQ(expect_iflagV,iflag_);
      }
      ASSERT_EQ(expectedRankV,rankVW);
      // check wether this worked out
      lidx_t ldaV;
      ST* V_vp=NULL;
      SUBR(mvec_extract_view)(V,&V_vp,&ldaV,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int nvec_V;
      SUBR(mvec_num_vectors)(V,&nvec_V,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(VTest::nvec_,nvec_V);

      lidx_t nloc_V;
      SUBR(mvec_my_length)(V,&nloc_V,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_,nloc_V);
      
      SUBR(mvec_from_device)(V,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),VTest::ColsAreNormalized(V_vp,nloc_,ldaV,stride_,mpi_comm_),tolV);
      ASSERT_NEAR(mt::one(),VTest::ColsAreOrthogonal(V_vp,nloc_,ldaV,stride_,mpi_comm_),tolV);
      
      int nsteps=2;

      // now orthogonalize W against V. The result should be such that Q*R1=W-V*R2, Q'*Q=I,V'*Q=0
      bool usedSVQB=false;
      rankVW=-42;
      SUBR(orthog)(V,Q,NULL,R1,R2,nsteps,&rankVW,&iflag_);
      if (iflag_!=+2)
      {
        ASSERT_EQ(expect_iflagVW,iflag_);
      } 
      else
      {
        // +2 means mvec_QR returned -99 (not implemented), and SVQB was used instead
        // so that the matrix R1 is not an upper triangular matrix but W is still made
        // orthogonal to V
        usedSVQB=true;
      }
      ASSERT_EQ(expectedRankVW,rankVW);
      
      // check orthonormality of Q
      lidx_t ldaQ;
      ST* Q_vp=NULL;
      SUBR(mvec_extract_view)(Q,&Q_vp,&ldaQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int nvec_Q;
      SUBR(mvec_num_vectors)(Q,&nvec_Q,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(WTest::nvec_,nvec_Q);

      lidx_t nloc_Q;
      SUBR(mvec_my_length)(V,&nloc_Q,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_,nloc_Q);
      
      SUBR(mvec_from_device)(Q,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),WTest::ColsAreNormalized(Q_vp,nloc_,ldaQ,stride_,mpi_comm_),tolW);
      ASSERT_NEAR(mt::one(),WTest::ColsAreOrthogonal(Q_vp,nloc_,ldaQ,stride_,mpi_comm_),tolW);

      if (usedSVQB)
      {
        //TODO - probably in this case we could ensure that the relation
        //       Q=(W-V*R2)*B [with R1=B] holds, but currently orthog does
        //       not promise this so we don't check for it.
      }
      else
      {
        // check the decomposition: Q*R1 = W - V*R2 (compute W2=Q*R1+V*R2-W and compare with 0)
        SUBR(mvec_times_sdMat)(st::one(),Q,R1,st::zero(),W2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_times_sdMat)(st::one(),V,R2,st::one(),W2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_add_mvec)(-st::one(),W,st::one(),W2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_from_device)(W2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),ArrayEqual(W2_vp_,nloc_,k_,ldaW2_,stride_,st::zero(),vflag_),tolW);
      }
  }

};

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
      doOrthogTests(V_, W_, Q_, R0_, R1_, R2_, 
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
      doOrthogTests(V_, W_, Q_, R0_, R1_, R2_, 
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
      
      lidx_t nR=_M_+_K_;
      TYPE(sdMat_ptr) R=NULL;
      SUBR(sdMat_create)(&R,nR,nR,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(R0_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R0_=NULL;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=NULL;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=NULL;
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
      doOrthogTests(V_, W_, Q_, R0_, R1_, R2_, 
        0, 0, m_, m_+k_);

      // check the location of the resulting R parts
      _ST_ *R_vp;
      lidx_t ldaR;

      SUBR(sdMat_extract_view)(R,&R_vp,&ldaR,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

#ifdef PHIST_SDMATS_ROW_MAJOR
#warning "test not implemented for row-major sdMats"
return;
#endif

        _MT_ errR2=mt::one();
        for (lidx_t i=0; i<m_; i++)
        {
          for (lidx_t j=0; j<k_; j++)
          {
            lidx_t jj = m_+j;
            _ST_ a = R_vp[jj*ldaR+i];
            _ST_ b = R2_vp_[j*ldaR2_+i];
            _MT_ denom=st::abs(a+b);
            if (denom==mt::zero()) denom=mt::one();
            errR2+= st::abs(a-b)/denom;
          }
        }
        _MT_ errR1=mt::one();
        for (lidx_t i=0; i<k_; i++)
        {
          for (lidx_t j=0; j<k_; j++)
          {
            lidx_t ii = m_+i;
            lidx_t jj = m_+j;
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
      R0_=NULL;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=NULL;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=NULL;
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
      
      lidx_t nR=_M_+_K_;
      TYPE(sdMat_ptr) R=NULL;
      SUBR(sdMat_create)(&R,nR,nR,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(R0_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R0_=NULL;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=NULL;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=NULL;
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
      doOrthogTests(V_, W_, Q_, R0_, R1_, R2_, 
        m_>1? 1:0, 1, 1, m_);
      
      // check the location of the resulting R parts
      _ST_ *R_vp;
      lidx_t ldaR;

      SUBR(sdMat_extract_view)(R,&R_vp,&ldaR,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);

#ifdef PHIST_SDMATS_ROW_MAJOR
#warning "test not implemented for row-major sdMats"
        return;
#endif

        _MT_ errR2=mt::one();
        for (lidx_t i=0; i<m_; i++)
        {
          for (lidx_t j=0; j<k_; j++)
          {
            lidx_t jj = m_+j;
            _ST_ a = R_vp[jj*ldaR+i];
            _ST_ b = R2_vp_[j*ldaR2_+i];
            _MT_ denom=st::abs(a+b);
            if (denom==mt::zero()) denom=mt::one();
            errR2+= st::abs(a-b)/denom;
          }
        }
        _MT_ errR1=mt::one();
        for (lidx_t i=0; i<k_; i++)
        {
          for (lidx_t j=0; j<k_; j++)
          {
            lidx_t ii = m_+i;
            lidx_t jj = m_+j;
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
      R0_=NULL;
      SUBR(sdMat_delete)(R1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R1_=NULL;
      SUBR(sdMat_delete)(R2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      R2_=NULL;
      SUBR(sdMat_delete)(R,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check if we can orthogonalize W against V where V and W are views into Z=[V W].
  TEST_F(CLASSNAME, random_vectors_viewing_same_block)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      TYPE(mvec_ptr) V_big=NULL, V=NULL, W=NULL;
      int ncols = m_+k_+13;// some extra padding to make the test more interesting
      PHISTTEST_MVEC_CREATE(&V_big,map_,ncols,&iflag_);
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
      doOrthogTests(V, W, Q_, R0_, R1_, R2_, 
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
  TEST_F(CLASSNAME, one_vectors_viewing_same_block)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      TYPE(mvec_ptr) V_big=NULL, V=NULL, W=NULL;
      int ncols = m_+k_+13;// some extra padding to make the test more interesting
      PHISTTEST_MVEC_CREATE(&V_big,map_,ncols,&iflag_);
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
      doOrthogTests(V, W, Q_, R0_, R1_, R2_, 
        m_>1? 1:0, 1, 1, m_);

      SUBR(mvec_delete)(V,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(W,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(V_big,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
