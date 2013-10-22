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
    SUBR(sdMat_random)(mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
    SUBR(sdMat_put_value)(mat2_,(ST)42.0,&ierr_);
      ASSERT_EQ(0,ierr_);
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
      SUBR(sdMat_get_nrows)(mat1_,&n,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(nrows_, n); 
      SUBR(sdMat_get_ncols)(mat1_,&n,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(ncols_, n); 
      }
    }

  // view certain rows and columns of a serial dense matrix,
  // manipulate them and check the original matrix has changed.
  TEST_F(CLASSNAME, view_block)
    {
    if (typeImplemented_)
      {
      int stride=1;
      int imin=std::min(2,nrows_-1);
      int imax=std::min(4,nrows_-1);
      int jmin=std::min(1,ncols_-1);
      int jmax=std::min(3,ncols_-1);
      TYPE(sdMat_ptr) m1_view=NULL;
      SUBR(sdMat_view_block)(mat1_,&m1_view,imin,imax,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      _ST_* val_ptr;
      int lda;
      SUBR(sdMat_extract_view)(m1_view,&val_ptr,&lda,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(lda,m_lda_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+imin+jmin*lda,val_ptr,imax-imin+1,jmax-jmin+1,lda,stride));

      SUBR(sdMat_print)(mat1_,&ierr_);
      SUBR(sdMat_print)(m1_view,&ierr_);

      // set all the viewed entries to a certain value and check that the original vector is
      // changed.
      _ST_ val = random_number();
      SUBR(sdMat_put_value)(m1_view,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      SUBR(sdMat_print)(m1_view,&ierr_);
      
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_+imin+jmin*lda,imax-imin+1,jmax-jmin+1,lda,stride,val));
      }

    }

