#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public KernelTestWithSdMats<_ST_,_NV_,_NV_> 
  {

public:


  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::SetUp();
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::TearDown();
    }

void BuildTestCase1()
  {
  int ilower;
  phist_map_get_ilower(map_,&ilower,&ierr_);
  ASSERT_EQ(0,ierr_);
  for (int j=0;j<nvec_;j++)
    for (int i=0; i<stride_*nloc_; i+=stride_)
      {
      MT ij = (MT)((i+ilower+1)*(j+1));
      vec1_vp_[j*lda_+i] = ij*st::one() + ij*st::I();
      vec2_vp_[j*lda_+i] = st::one()/(ij - ij*st::I());
      }
  // result of vec1'*vec2
  for (int j=0; j<ncols_; j++)
    for (int i=0; i<nrows_; i++)
    {
    mat2_vp_[j*m_lda_+i] = (MT)((i+1)*nglob_)*st::one()/(MT)(j+1);
    }
  }

// this is just for debugging
void PrintTestCase()
  {
  std::cout << std::setw(8) << std::setprecision(8);
  std::cout << "i="<<st::I()<<std::endl;
  std::cout << "A="<<std::endl;
  for (int i=0; i<stride_*nloc_; i+=stride_)
    {
    for (int j=0;j<nvec_;j++)
      {
      std::cout << vec1_vp_[j*lda_+i] << "  ";
      }
    std::cout << std::endl;
    }
  std::cout << "B="<<std::endl;
  for (int i=0; i<stride_*nloc_; i+=stride_)
    {
    for (int j=0;j<nvec_;j++)
      {
      std::cout << vec2_vp_[j*lda_+i] << "  ";
      }
    std::cout << std::endl;
    }
  std::cout << "C=A'B: "<<std::endl;
  for (int i=0; i<nrows_; i++)
    {
    for (int j=0; j<ncols_; j++)
      {
      std::cout << mat2_vp_[j*m_lda_+i] <<"  ";
      }
    std::cout << std::endl;
    }
  }
};

  // vec1'*vec2 gives a real-valued matrix as defined in mat2
  TEST_F(CLASSNAME, with_real_result)
    {
    if (typeImplemented_)
      {
      BuildTestCase1();
//      PrintTestCase();
      _ST_ alpha=st::one(); 
      _ST_ beta=st::zero();
      SUBR(mvecT_times_mvec)(alpha,vec1_,vec2_,beta,mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,1));

      alpha=(_MT_)0.5*st::one();
      beta=(_MT_)0.5*st::one();
      SUBR(mvecT_times_mvec)(alpha,vec1_,vec2_,beta,mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,1));

      alpha=(_MT_)0.3*st::one(); 
      beta=(_MT_)0.7*st::one();
      SUBR(mvecT_times_mvec)(alpha,vec1_,vec2_,beta,mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,1));
      }
    }
    

