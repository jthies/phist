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
  
  //! M is n x n, M=S*T with T (block-) upper triangular
  TYPE(sdMat_ptr) M_,S_,T_;

  std::complex<MT> ev_[_N_];  
  
  /*! Set up routine.
   */
  virtual void SetUp()
    {
    MTest::SetUp();
    
    // fill all matrices with random numbers and set some pointers
    M_=mat1_;
    T_=mat1_; // SchurDecomp works in-place
    S_=mat2_;
    }

  /*! Clean up.
   */
  virtual void TearDown()
    {
    KernelTestWithType<_ST_>::TearDown();
    }

  void DoSchurDecomp(int nselect,int nsort,eigSort_t which)
    {

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
    }
};

  // check ones(n,m)'*ones(n,m)=n*ones(m,m)
  TEST_F(CLASSNAME, do_tests) 
    {
    if (typeImplemented_)
      {
      eigSort_t sort[4];
      sort[0]=LM;
      sort[1]=SM;
      sort[2]=LR;
      sort[3]=SR;
      
      int nselect[4];
      int nsort[4];
      // a few test cases
      nselect[0]=1; nsort[0]=1;
      nselect[1]=_N_; nsort[1]=_N_;
      nselect[2]=std::min(5,_N_); nsort[2]=0;
      nselect[3]=std::min(7,_N_); nsort[3]=3;
      for (int s=0;s<4;s++)
        {
        for (int c=0;c<4;c++)
          {
          PHIST_OUT(PHIST_INFO,"==================================================");
          PHIST_OUT(PHIST_INFO,"nselect %d, nsort %d, order=%s",nselect[c],nsort[c],
                sort[s]==LM?"LM":sort[s]==SM?"SM":sort[s]==LR?"LR":sort[s]==SR?"SR":"??");
          PHIST_OUT(PHIST_INFO,"==================================================");
          DoSchurDecomp(nselect[c],nsort[c],sort[s]);
          // check the T matrix is upper triangular (with 2x2 blocks for complex eigs in the 
          // real case)
          PHIST_DEB("data types: ST %c, MT %c, CT %c",
                st::type_char(), mt::type_char(), ct::type_char());
#ifdef _IS_COMPLEX_
          for (int j=0;j<_N_;j++)
            {
            for (int i=j+1;i<_N_;i++)
              {
              ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[j*m_lda_+i]));
              }//i
            }//j
#else
          for (int j=0;j<_N_;j++)
            {
            PHIST_DEB("j=%d, ev[j]=%8.4f %+8.4fi",j,ct::real(ev_[j]),ct::imag(ev_[j]));
            if (ct::imag(ev_[j])<=mt::zero())
              {
              ASSERT_REAL_EQ(st::zero(),st::abs(mat1_vp_[j*m_lda_+j+1]));
              }
            for (int i=j+2;i<_N_;i++)
              {
              ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[j*m_lda_+i]));
              }//i
            }//j
#endif
          }//c(ase)
        }//s(ort-type)
      }// type implemented?
    }
