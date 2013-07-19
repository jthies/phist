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

// note: the numerical functions of sdMat are not tested up to now. This makes sense as
// in all three kernel libs we want to support (Epetra, Tpetra, Ghost) an sdMat is just
// an mvec with a local map.
  TEST_F(CLASSNAME, get_attributes)
    {
    if (typeImplemented_)
      {
      int n;
      _SUBR_(sdMat_get_nrows)(mat1_,&n,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(nrows_, n); 
      _SUBR_(sdMat_get_ncols)(mat1_,&n,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(ncols_, n); 
      }
    }

