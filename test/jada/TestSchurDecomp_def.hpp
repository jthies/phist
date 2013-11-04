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
    SUBR(sdMat_random)(mat1_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_random)(mat2_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
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
    SUBR(SchurDecomp)(mat1_vp_,m_lda_,mat2_vp_,m_lda_,n_,nselect,nsort,which,ev_,&this->ierr_);
    ASSERT_EQ(0,ierr_);
    }
};

  // check ones(n,m)'*ones(n,m)=n*ones(m,m)
  TEST_F(CLASSNAME, do_tests) 
    {
    if (typeImplemented_)
      {
      eigSort_t sort[4];
      sort[0]=SR;
      sort[1]=LR;
      sort[2]=SM;
      sort[3]=LM;
      
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
          DoSchurDecomp(nselect[c],nsort[c],sort[s]);
          // check the T matrix is upper triangular (with 2x2 blocks for complex eigs in the 
          // real case)
#ifdef _IS_COMPLEX_
          for (int i=0;i<_N_;i++)
            {
            for (int j=i+1;j<_N_;j++)
              {
              ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[j*m_lda_+i]));
              }//i
            }//j
#else
          for (int j=0;j<_N_;j++)
            {
            if (ct::imag(ev_[j])==mt::zero())
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
