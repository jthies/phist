#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_N_,_N_>
  {

public:

  typedef KernelTestWithSdMats<_ST_,_N_,_N_> MTest;

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_N_;

  // different test scenarios
  int nselect_[4];
  int nsort_[4];

  std::complex<MT> ev_[_N_];
  
  /*! Set up routine.
   */
  virtual void SetUp()
    {
    MTest::SetUp();
      
    // a few test cases
    nselect_[0]=1; nsort_[0]=1;
    nselect_[1]=_N_; nsort_[1]=_N_;
    nselect_[2]=std::min(5,_N_); nsort_[2]=0;
    nselect_[3]=std::min(7,_N_); nsort_[3]=3;
    }

  /*! Clean up.
   */
  virtual void TearDown()
    {
    KernelTestWithType<_ST_>::TearDown();
    }

  void DoSchurDecompTest(eigSort_t which)
    {
    if (!typeImplemented_) return;
    for (int c=0;c<4;c++) 
      {
    int nselect = nselect_[c];
    int nsort = nsort_[c];
    PHIST_OUT(PHIST_INFO,"==================================================");
    PHIST_OUT(PHIST_INFO,"CASE nselect %d, nsort %d",nselect,nsort);
    PHIST_OUT(PHIST_INFO,"==================================================");
    
    SUBR(sdMat_random)(mat1_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    
    SUBR(sdMat_random)(mat2_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);

    PHIST_DEB("input matrix to Schur-decomp:");
    SUBR(sdMat_print)(mat1_,&ierr_);
    SUBR(SchurDecomp)(mat1_vp_,m_lda_,mat2_vp_,m_lda_,n_,nselect,nsort,which,ev_,&this->ierr_);
    PHIST_DEB("resulting T:");
    SUBR(mvec_print)(mat1_,&ierr_);
    ASSERT_EQ(0,ierr_);
    
    PHIST_DEB("eigenvalue array");
    for (int i=0;i<n_;i++)
      {
      // test the traits class on the way:
      ASSERT_REAL_EQ(ct::abs(ev_[i]),std::abs(ev_[i]));
      ASSERT_REAL_EQ(ct::real(ev_[i]),std::real(ev_[i]));
      ASSERT_REAL_EQ(ct::imag(ev_[i]),std::imag(ev_[i]));
      PHIST_DEB("%8.4f%+8.4fi\tabs=%8.4f",ct::real(ev_[i]),ct::imag(ev_[i]),ct::abs(ev_[i]));
      }
    
    // check that the eigenvalues on the diagonal of T have the same ordering as those in 
    // ev_
    for (int i=0;i<n_;i++)
      {
      ASSERT_REAL_EQ(ct::real(ev_[i]), st::real(mat1_vp_[i*m_lda_+i]));
#ifdef _IS_COMPLEX_
      ASSERT_REAL_EQ(ct::imag(ev_[i]), st::imag(mat1_vp_[i*m_lda_+i]));
#endif      
      }
    // check that the first nsort eigenvalues are sorted correctly
    for (int i=1;i<nsort;i++)
      {
      if (which==LM)
        {
        ASSERT_TRUE(ct::abs(ev_[i])<=ct::abs(ev_[i-1]));
        }
      else if (which==SM)
        {
        ASSERT_TRUE(ct::abs(ev_[i])>=ct::abs(ev_[i-1]));
        }      
      else if (which==LR)
        {
        ASSERT_TRUE(ct::real(ev_[i])<=ct::real(ev_[i-1]));
        }
      else if (which==SR)
        {
        ASSERT_TRUE(ct::real(ev_[i])>=ct::real(ev_[i-1]));
        }
      }
    // check that the first nselect are the largest/smallest etc globally
    MT val;
    if (which==LM||which==SM) val=ct::abs(ev_[0]);
    else if (which==LR||which==SR) val=ct::real(ev_[0]);
    for (int i=1;i<nselect;i++)
      {
      if (which==LM)
        {
        val=std::min(ct::abs(ev_[i]),val);
        }
      else if (which==SM)
        {
        val=std::max(ct::abs(ev_[i]),val);
        }      
      else if (which==LR)
        {
        val=std::min(ct::real(ev_[i]),val);
        }
      else if (which==SR)
        {
        val=std::max(ct::real(ev_[i]),val);
        }
      }
    MT err=(MT)0.0;
    for (int i=nselect;i<n_;i++)
      {
      if (which==LM)
        {
        // make sure val >= all others in abs value
        err=std::max(err,val-ct::abs(ev_[i]));
        }
      else if (which==SM)
        {
        // make sure val<=all others in abs value
        err=std::max(err,ct::abs(ev_[i])-val);
        }
      if (which==LR)
        {
        err=std::max(err,val-ct::real(ev_[i]));
        }
      else if (which==SR)
        {
        err=std::max(err,ct::real(ev_[i])-val);
        }
      }
    // make sure err>0
    ASSERT_REAL_EQ(mt::zero(),std::min(mt::zero(),err));
    
      // check the T matrix is upper triangular (with 2x2 blocks for complex eigs in the 
      // real case)
#ifdef _IS_COMPLEX_
      for (int j=0;j<_N_;j++)
        {
        for (int i=j+1;i<_N_;i++)
          {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+i]));
          }//i
        }//j
#else
      for (int j=0;j<_N_;j++)
        {
        PHIST_DEB("j=%d, ev[j]=%8.4f %+8.4fi",j,ct::real(ev_[j]),ct::imag(ev_[j]));
        if (ct::imag(ev_[j])<=mt::zero() && j+1<_N_)
          {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+j+1]));
          }
        for (int i=j+2;i<_N_;i++)
          {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+i]));
          }//i
        }//j
#endif
      }//cases
    }//DoSchurDecompTest
  };

  TEST_F(CLASSNAME, rand_LM) 
    {
    DoSchurDecompTest(LM);
    }

  TEST_F(CLASSNAME, rand_SM) 
    {
    DoSchurDecompTest(SM);
    }

  TEST_F(CLASSNAME, rand_LR) 
    {
    DoSchurDecompTest(LR);
    }

  TEST_F(CLASSNAME, rand_SR) 
    {
    DoSchurDecompTest(SR);
    }

