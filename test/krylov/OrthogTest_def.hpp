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
  _TYPE_(mvec_ptr) V_, W_;

  //! R1 is m x m, R2 is k x k
  _TYPE_(sdMat_ptr) R1_, R2_;
  
  _ST_ *V_vp_,*W_vp_,*R1_vp_,*R2_vp_;
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  int ldaV_,ldaW_,ldaR1_,ldaR2_,stride_;

  
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
      _SUBR_(mvec_create)(&V_,this->map_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_extract_view)(V_,&V_vp_,&ldaV_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_create)(&W_,this->map_,this->k_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_extract_view)(W_,&W_vp_,&ldaW_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      // create matrices R1, R2 and matrix views for setting/checking entries
      _SUBR_(sdMat_create)(&R1_,this->m_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(sdMat_extract_view)(R1_,&R1_vp_,&this->ldaR1_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(sdMat_create)(&R2_,this->k_,this->k_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_extract_view)(R2_,&R2_vp_,&this->ldaR2_,&this->ierr_);
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
      _SUBR_(mvec_delete)(V_,&ierr_);
      _SUBR_(mvec_delete)(W_,&ierr_);
      _SUBR_(sdMat_delete)(R1_,&ierr_);
      _SUBR_(sdMat_delete)(R2_,&ierr_);
      }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
    }

};

  // check if vectors are normalized correctly after QR factorization
  TEST_F(CLASSNAME, do_some_test) 
    {
    if (typeImplemented_)
      {
      // fill V and W with random numbers
      _SUBR_(mvec_random)(V_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // orthogonalize the m columns of V
      _SUBR_(mvec_QR)(V_,R1_,&ierr_);
      ASSERT_EQ(0,ierr_);

      std::cout << "V="<<std::endl;
      for (int i=0;i<n_;i++)
        {
        std::cout << i << "\t";
        for (int j=0;j<m_;j++)
          {
          std::cout <<  V_vp_[j*ldaV_+i]<<"  ";
          }
        std::cout << std::endl;
        }
      // check wether this worked out
      ASSERT_REAL_EQ((_MT_)1.0,VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_));
      ASSERT_REAL_EQ((_MT_)1.0,VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_));
      // and so on.
      }
    }


