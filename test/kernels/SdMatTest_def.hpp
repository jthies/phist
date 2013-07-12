#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_> 
  {

public:


  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithSdMats<ST,_NROWS_,_NCOLS_>::SetUp();
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithSdMats<ST,_NROWS_,_NCOLS_>::TearDown();
    }

};

  TEST_F(CLASSNAME, put_value) 
    {
    if (typeImplemented_)
      {
      ST val = 42.0 + 3.0*mt::I();
      _SUBR_(sdMat_put_value)(mat1_,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_,nrows_,ncols_,m_lda_,1,val));
      }
    }
    

