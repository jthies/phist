#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_> 
  {

public:


  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::SetUp();
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
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
  gidx_t ilower;
  phist_map_get_ilower(map_,&ilower,&iflag_);
  ASSERT_EQ(0,iflag_);
  for (int j=0;j<nvec_;j++)
    for (int i=0; i<stride_*nloc_; i+=stride_)
      {
      MT ij = (MT)((i+ilower+1)*(j+1));
#ifdef PHIST_MVECS_ROW_MAJOR
      vec1_vp_[j+i*lda_] = ij*st::one() + ij*st::cmplx_I();
      vec2_vp_[j+i*lda_] = st::one()/(ij - ij*st::cmplx_I());
#else
      vec1_vp_[j*lda_+i] = ij*st::one() + ij*st::cmplx_I();
      vec2_vp_[j*lda_+i] = st::one()/(ij - ij*st::cmplx_I());
#endif
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
  std::cout << "i="<<st::cmplx_I()<<std::endl;
  std::cout << "A="<<std::endl;
  for (int i=0; i<stride_*nloc_; i+=stride_)
    {
    for (int j=0;j<nvec_;j++)
      {
#ifdef PHIST_MVECS_ROW_MAJOR
      std::cout << vec1_vp_[j+i*lda_] << "  ";
#else
      std::cout << vec1_vp_[j*lda_+i] << "  ";
#endif
      }
    std::cout << std::endl;
    }
  std::cout << "B="<<std::endl;
  for (int i=0; i<stride_*nloc_; i+=stride_)
    {
    for (int j=0;j<nvec_;j++)
      {
#ifdef PHIST_MVECS_ROW_MAJOR
      std::cout << vec2_vp_[j+i*lda_] << "  ";
#else
      std::cout << vec2_vp_[j*lda_+i] << "  ";
#endif
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
    if (typeImplemented_ && !problemTooSmall_)
      {
      BuildTestCase1();
//      PrintTestCase();
      _ST_ alpha=st::one(); 
      _ST_ beta=st::zero();
      SUBR(mvecT_times_mvec)(alpha,vec1_,vec2_,beta,mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat1_,mat2_),100*mt::eps());

      alpha=(_MT_)0.5*st::one();
      beta=(_MT_)0.5*st::one();
      SUBR(mvecT_times_mvec)(alpha,vec1_,vec2_,beta,mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat1_,mat2_),100*mt::eps());

      alpha=(_MT_)0.3*st::one(); 
      beta=(_MT_)0.7*st::one();
      SUBR(mvecT_times_mvec)(alpha,vec1_,vec2_,beta,mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat1_,mat2_), 100*mt::eps());
      }
    }
    

