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
    if( typeImplemented_ )
    {
      SUBR(sdMat_random)(mat1_,&ierr_);
        ASSERT_EQ(0,ierr_);
      SUBR(sdMat_put_value)(mat2_,(ST)42.0,&ierr_);
        ASSERT_EQ(0,ierr_);
      SUBR(sdMat_random)(mat3_,&ierr_);
        ASSERT_EQ(0,ierr_);
    }
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

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(mat1_,&ierr_);
      SUBR(sdMat_print)(m1_view,&ierr_);
#endif
      
      _ST_* val_ptr;
      int lda;
      SUBR(sdMat_extract_view)(m1_view,&val_ptr,&lda,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(lda,m_lda_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+imin+jmin*lda,val_ptr,imax-imin+1,jmax-jmin+1,lda,stride));

      // set all the viewed entries to a certain value and check that the original vector is
      // changed.
      _ST_ val = random_number();
      SUBR(sdMat_put_value)(m1_view,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(m1_view,&ierr_);
#endif
      
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_+imin+jmin*lda,imax-imin+1,jmax-jmin+1,lda,stride,val));
      }
    }

  // create a view of a view and check if this behaves as the user would suspect
  // (not knowing wether a sdmat is actually a view or not!)
  TEST_F(CLASSNAME, nested_view_block)
  {
    if (typeImplemented_)
    {
      // first set some data of the whole array
      _ST_ outer_val = st::rand();
      SUBR(sdMat_put_value)(mat1_,outer_val,&ierr_);
      ASSERT_EQ(0,ierr_);

      // now create a view
      int imin=std::min(3,nrows_-1);
      int imax=std::min(7,nrows_-1);
      int jmin=std::min(2,ncols_-1);
      int jmax=std::min(5,ncols_-1);
      TYPE(sdMat_ptr) view = NULL;
      SUBR(sdMat_view_block)(mat1_,&view,imin,imax,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);

      // set the data in the view to some other value
      _ST_ view_val = st::rand();
      SUBR(sdMat_put_value)(view,view_val,&ierr_);
      ASSERT_EQ(0,ierr_);

      // view part of view
      int nrows_view = imax-imin+1;
      int ncols_view = jmax-jmin+1;
      ASSERT_EQ(0,ierr_);
      int imin2=std::min(2,nrows_view-1);
      int imax2=std::min(4,nrows_view-1);
      int jmin2=std::min(1,ncols_view-1);
      int jmax2=std::min(3,ncols_view-1);
      TYPE(sdMat_ptr) view2 = NULL;
      SUBR(sdMat_view_block)(view, &view2, imin2,imax2,jmin2, jmax2, &ierr_);
      ASSERT_EQ(0,ierr_);

      // set data in the inner view to yet another value
      _ST_ inner_val = st::rand();
      SUBR(sdMat_put_value)(view2, inner_val, &ierr_);
      ASSERT_EQ(0,ierr_);

      // now use the raw data to verify results
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i < imin || i > imax || j < jmin || j > jmax )
          {
            ASSERT_REAL_EQ(mt::zero(), st::abs(mat1_vp_[j*m_lda_+i]-outer_val));
          }
          else if( i < imin+imin2 || i > imin+imax2 || j < jmin+jmin2 || j > jmin+jmax2 )
          {
            ASSERT_REAL_EQ(mt::zero(), st::abs(mat1_vp_[j*m_lda_+i]-view_val));
          }
          else
          {
            ASSERT_REAL_EQ(mt::zero(), st::abs(mat1_vp_[j*m_lda_+i]-inner_val));
          }
        }
      }
    }
  }


#ifdef IS_COMPLEX
  TEST_F(CLASSNAME, complex_random)
  {
    if(typeImplemented_)
    {
      // check that the imaginary part is not zero everywhere!
      MT maxAbsIm = mt::zero();
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols_; j++)
          maxAbsIm = std::max(std::abs(st::imag(mat1_vp_[j*m_lda_+i])), maxAbsIm);
      ASSERT_TRUE(maxAbsIm != mt::zero());
    }
  }
#endif

  // copy certain rows and columns of a serial dense matrix,
  // manipulate them and check the original matrix has not changed.
  TEST_F(CLASSNAME, get_set_block)
    {
    if (typeImplemented_)
      {
      int stride=1;
      int imin=std::min(2,nrows_-1);
      int imax=std::min(4,nrows_-1);
      int jmin=std::min(1,ncols_-1);
      int jmax=std::min(3,ncols_-1);
      TYPE(sdMat_ptr) m1_copy=NULL;
      
      // set m2=m1 to check later
      SUBR(sdMat_add_sdMat)(st::one(),mat1_,st::zero(),mat2_,&ierr_);
      
      SUBR(sdMat_create)(&m1_copy,imax-imin+1,jmax-jmin+1,comm_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(sdMat_get_block)(mat1_,m1_copy,imin,imax,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      _ST_* val_ptr;
      int lda;
      SUBR(sdMat_extract_view)(m1_copy,&val_ptr,&lda,&ierr_);
      ASSERT_EQ(0,ierr_);
      //note: can't use ArraysEqual here because we (may) have different strides      
      _MT_ maxval=mt::zero();
      for (int i=imin;i<=imax;i++)
        {
        for (int j=jmin;j<jmax;j++)
          {
          ST val1=mat1_vp_[j*m_lda_+i];
          ST val2=val_ptr[(j-jmin)*lda+(i-imin)];
          MT m = st::abs(val1-val2);
          MT p = st::abs(val1+val2);
          if (p==mt::zero()) p=mt::one();
          maxval=std::max(m/p,maxval);
          }
        }
      ASSERT_REAL_EQ(mt::one()+maxval,mt::one());

      // set all the copied entries to a certain value and check that the original vector 
      // is not changed.
      _ST_ val = random_number();
      SUBR(sdMat_put_value)(m1_copy,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));

      // copy back in the changed block
      SUBR(sdMat_set_block)(mat1_,m1_copy,imin,imax,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // check that the corresponding entries have changed
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_+imin+jmin*m_lda_,imax-imin+1,jmax-jmin+1,m_lda_,stride,val));
      }
    }

  TEST_F(CLASSNAME, sdMat_times_sdMat)
  {
    if (typeImplemented_ && nrows_ == ncols_ )
    {
      SUBR(sdMat_times_sdMat)(st::one(),mat1_,mat3_,st::one(),mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      PHIST_DEB("random*random+42");
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#endif

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols_; j++)
          for(int k = 0; k < ncols_; k++)
            mat2_vp_[j*m_lda_+i] -= mat1_vp_[k*m_lda_+i]*mat3_vp_[j*m_lda_+k];
      // check result
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,(ST)42.0),10*mt::eps());

    }
  }

  TEST_F(CLASSNAME, sdMatT_times_sdMat)
  {
    if (typeImplemented_ && nrows_ == ncols_ )
    {
#ifdef IS_COMPLEX
      // force some complex non-zero value even if complex_random fails
      mat1_vp_[0] = std::complex<MT>( st::imag(mat1_vp_[0]), st::real(mat1_vp_[0]) );
#endif
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat3_,st::one(),mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      PHIST_DEB("random'*random+42");
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#endif

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols_; j++)
          for(int k = 0; k < ncols_; k++)
            mat2_vp_[j*m_lda_+i] -= st::conj(mat1_vp_[i*m_lda_+k])*mat3_vp_[j*m_lda_+k];
      // check result
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,(ST)42.0),10*mt::eps());

    }
  }
