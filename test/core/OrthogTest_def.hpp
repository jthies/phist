#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithType< _ST_ >,
                 public KernelTestWithMap<_N_>
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
  int ldaV_,ldaW_,ldaW2_,ldaQ_,ldaR0_,ldaR1_,ldaR2_,stride_;

  
  //NOTE: we assume stride_=1 here to make the loops simpler,
  //      it seems reasonable to me to do that because mvecs 
  //      are generally stored in column-major order.

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithType<_ST_>::SetUp();
    KernelTestWithMap<_N_>::SetUp();
    if (this->typeImplemented_)
      {
      // create vectors V, W and vector views for setting/checking entries
      SUBR(mvec_create)(&V_,this->map_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_extract_view)(V_,&V_vp_,&ldaV_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_create)(&W_,this->map_,this->k_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_extract_view)(W_,&W_vp_,&ldaW_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_create)(&Q_,this->map_,this->k_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_extract_view)(Q_,&Q_vp_,&ldaQ_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_create)(&W2_,this->map_,this->k_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_extract_view)(W2_,&W2_vp_,&ldaW2_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      // create matrices R0,R1, R2 and matrix views for setting/checking entries
      SUBR(sdMat_create)(&R0_,this->m_,this->m_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(sdMat_extract_view)(R0_,&R0_vp_,&this->ldaR0_,&this->ierr_);
      SUBR(sdMat_create)(&R1_,this->k_,this->k_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(sdMat_extract_view)(R1_,&R1_vp_,&this->ldaR1_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(sdMat_create)(&R2_,this->m_,this->k_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(sdMat_extract_view)(R2_,&R2_vp_,&this->ldaR2_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      }
    stride_=1;
    }

  /*! Clean up.
   */
  virtual void TearDown()
    {
    if (this->typeImplemented_)
      {
      SUBR(mvec_delete)(V_,&ierr_);
      SUBR(mvec_delete)(W_,&ierr_);
      SUBR(mvec_delete)(W2_,&ierr_);
      SUBR(mvec_delete)(Q_,&ierr_);
      SUBR(sdMat_delete)(R0_,&ierr_);
      SUBR(sdMat_delete)(R1_,&ierr_);
      SUBR(sdMat_delete)(R2_,&ierr_);
      }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
    }

};

  // check if vectors are normalized correctly after QR factorization
  TEST_F(CLASSNAME, test_with_random_vectors) 
    {
    if (typeImplemented_)
      {
      // fill V and W with random numbers
      SUBR(mvec_random)(V_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_random)(W_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // copy Q=W because orthog() works in-place
      SUBR(mvec_add_mvec)(st::one(),W_,st::zero(),Q_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // orthogonalize the m columns of V
      SUBR(mvec_QR)(V_,R0_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // check wether this worked out
      ASSERT_REAL_EQ(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_));
      ASSERT_REAL_EQ(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_));
      
      int nsteps=2;

      // now orthogonalize W against V. The result should be such that Q*R1=W-V*R2, Q'*Q=I,V'*Q=0
      SUBR(orthog)(V_,Q_,R1_,R2_,nsteps,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // check orthonormality of Q
      ASSERT_REAL_EQ(mt::one(),WTest::ColsAreNormalized(Q_vp_,nloc_,ldaQ_,stride_,mpi_comm_));
      ASSERT_REAL_EQ(mt::one(),WTest::ColsAreOrthogonal(Q_vp_,nloc_,ldaQ_,stride_,mpi_comm_));
      
      // check the decomposition: Q*R2 = W - V*R1 (compute W2=Q*R2+V*R1-W and compare with 0)
      SUBR(mvec_times_sdMat)(st::one(),Q_,R1_,st::zero(),W2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_times_sdMat)(st::one(),V_,R2_,st::one(),W2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(-st::one(),W_,st::one(),W2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_NEAR(mt::one(),ArrayEqual(W2_vp_,nloc_,k_,ldaW2_,stride_,st::zero()),8*mt::eps());
      }
    }

