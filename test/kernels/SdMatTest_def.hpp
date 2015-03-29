#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_> 
  {

public:

  typedef KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_> MTest;

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    KernelTestWithSdMats<ST,_NROWS_,_NCOLS_>::SetUp();
    if( typeImplemented_ )
    {
      SUBR(sdMat_random)(mat1_,&iflag_);
        ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(mat2_,(ST)42.0,&iflag_);
        ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(mat3_,&iflag_);
        ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithSdMats<ST,_NROWS_,_NCOLS_>::TearDown();
    }

};

  TEST_F(CLASSNAME, get_attributes)
    {
    if (typeImplemented_)
      {
      int n;
      SUBR(sdMat_get_nrows)(mat1_,&n,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nrows_, n); 
      SUBR(sdMat_get_ncols)(mat1_,&n,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(ncols_, n); 
      }
    }

  TEST_F(CLASSNAME, upload_download)
  {
    // just tests that the upload and download functions return 0
    // and do not crash.
    if (typeImplemented_)
    {
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // X = 1*Y + 0*X = Y
  TEST_F(CLASSNAME, copy_by_axpy)
  {
    if (typeImplemented_)
    {
      ST alpha = st::one();
      ST beta  = st::zero();
      SUBR(sdMat_add_sdMat)(alpha,mat1_,beta,mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
    }

  }

  // X = 0*Y + a*X = a*X
  TEST_F(CLASSNAME, scale_by_axpy)
  {
    if (typeImplemented_)
    {
      int stride = 1;
      ST alpha = st::zero();
      ST beta  = st::prand();
      PHIST_OUT(9,"axpy, alpha=%f+%f i, beta=%f+%f i\n",st::real(alpha),
          st::imag(alpha),st::real(beta),st::imag(beta));
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      MTest::PrintSdMat(std::cerr,"before scale",mat2_vp_,m_lda_,stride,mpi_comm_);
#endif
      SUBR(sdMat_add_sdMat)(alpha,mat1_,beta,mat2_,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      MTest::PrintSdMat(std::cerr,"after scale",mat2_vp_,m_lda_,stride,mpi_comm_);
#endif
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),SdMatEqual(mat2_,beta*((ST)42.0)));
    }
  }

  // X = a*Y + b*X
  TEST_F(CLASSNAME, random_add)
  {
    if( typeImplemented_ )
    {
      ST beta = st::prand();
      ST alpha= st::prand();

      SUBR(sdMat_add_sdMat)(alpha,mat1_,beta,mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // calculate solution by hand
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          mat1_vp_[MIDX(i,j,m_lda_)] = alpha*mat1_vp_[MIDX(i,j,m_lda_)]+42.0*beta;
        }
      }
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
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
      SUBR(sdMat_view_block)(mat1_,&m1_view,imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      _ST_* val_ptr;
      lidx_t lda;
      SUBR(sdMat_extract_view)(m1_view,&val_ptr,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda,m_lda_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+imin+jmin*lda,val_ptr,imax-imin+1,jmax-jmin+1,lda,stride));

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(mat1_,&iflag_);
      SUBR(sdMat_print)(m1_view,&iflag_);
#endif
      
      SUBR(sdMat_extract_view)(m1_view,&val_ptr,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda,m_lda_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+MIDX(imin,jmin,lda),val_ptr,imax-imin+1,jmax-jmin+1,lda,stride,mflag_));

      // set all the viewed entries to a certain value and check that the original vector is
      // changed.
      _ST_ val = random_number();
      SUBR(sdMat_put_value)(m1_view,val,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+imin+jmin*lda,val_ptr,imax-imin+1,jmax-jmin+1,lda,stride));
      
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(m1_view,&iflag_);
#endif
      
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(&mat1_vp_[MIDX(imin,jmin,lda)],imax-imin+1,jmax-jmin+1,lda,stride,val,mflag_));

      SUBR(sdMat_delete)(m1_view,&iflag_);
      ASSERT_EQ(0,iflag_);
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
      SUBR(sdMat_put_value)(mat1_,outer_val,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // now create a view
      int imin=std::min(3,nrows_-1);
      int imax=std::min(7,nrows_-1);
      int jmin=std::min(2,ncols_-1);
      int jmax=std::min(5,ncols_-1);
      TYPE(sdMat_ptr) view = NULL;
      SUBR(sdMat_view_block)(mat1_,&view,imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      // set the data in the view to some other value
      _ST_ view_val = st::rand();
      SUBR(sdMat_put_value)(view,view_val,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // view part of view
      int nrows_view = imax-imin+1;
      int ncols_view = jmax-jmin+1;
      ASSERT_EQ(0,iflag_);
      int imin2=std::min(2,nrows_view-1);
      int imax2=std::min(4,nrows_view-1);
      int jmin2=std::min(1,ncols_view-1);
      int jmax2=std::min(3,ncols_view-1);
      TYPE(sdMat_ptr) view2 = NULL;
      SUBR(sdMat_view_block)(view, &view2, imin2,imax2,jmin2, jmax2, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      // set data in the inner view to yet another value
      _ST_ inner_val = st::rand();
      SUBR(sdMat_put_value)(view2, inner_val, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);

      // now use the raw data to verify results
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i < imin || i > imax || j < jmin || j > jmax )
          {
            ASSERT_REAL_EQ(mt::zero(), st::abs(mat1_vp_[MIDX(i,j,m_lda_)]-outer_val));
          }
          else if( i < imin+imin2 || i > imin+imax2 || j < jmin+jmin2 || j > jmin+jmax2 )
          {
            ASSERT_REAL_EQ(mt::zero(), st::abs(mat1_vp_[MIDX(i,j,m_lda_)]-view_val));
          }
          else
          {
            ASSERT_REAL_EQ(mt::zero(), st::abs(mat1_vp_[MIDX(i,j,m_lda_)]-inner_val));
          }
        }
      }

      SUBR(sdMat_delete)(view2,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(view,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


#ifdef IS_COMPLEX
  TEST_F(CLASSNAME, complex_random)
  {
    if(typeImplemented_)
    {
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check that the imaginary part is not zero everywhere!
      MT maxAbsIm = mt::zero();
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols_; j++)
          maxAbsIm = std::max(std::abs(st::imag(mat1_vp_[MIDX(i,j,m_lda_)])), maxAbsIm);
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
      SUBR(sdMat_add_sdMat)(st::one(),mat1_,st::zero(),mat2_,&iflag_);

      // did the assignment operation above work, M2=M1?
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
      
      SUBR(sdMat_create)(&m1_copy,imax-imin+1,jmax-jmin+1,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_get_block)(mat1_,m1_copy,imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);
      
#if PHIST_OUTLEV>=PHIST_DEBUG && !defined(PHIST_KERNEL_LIB_BUILTIN)
      PHIST_SOUT(PHIST_DEBUG,"i=[%" PRlidx ",%" PRlidx "], j=[%" PRlidx ",%" PRlidx "]\n",
        imin,imax,jmin,jmax);
      PHIST_SOUT(PHIST_DEBUG,"Original matrix:\n");
      SUBR(sdMat_print)(mat1_,&iflag_);
      PHIST_SOUT(PHIST_DEBUG,"This should be a copy of those columns:\n");
      SUBR(sdMat_print)(m1_copy,&iflag_);
#endif      
      _ST_* val_ptr;
      lidx_t lda;
      SUBR(sdMat_extract_view)(m1_copy,&val_ptr,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(m1_copy,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      //note: can't use ArraysEqual here because we (may) have different strides      
      _MT_ maxval=mt::zero();
      for (int i=imin;i<=imax;i++)
        {
        for (int j=jmin;j<jmax;j++)
          {
          ST val1=mat1_vp_[MIDX(i,j,m_lda_)];
          ST val2=val_ptr[MIDX(i-imin,j-jmin,lda)];
          MT mn = st::abs(val1-val2);
          MT pl = st::abs(val1+val2);
          if (pl==mt::zero()) pl=mt::one();
          maxval=std::max(mn/pl,maxval);
          }
        }
      ASSERT_REAL_EQ(mt::one()+maxval,mt::one());

      // set all the copied entries to a certain value and check that the original vector 
      // is not changed.
      _ST_ val = random_number();
      SUBR(sdMat_put_value)(m1_copy,val,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check that the pointers we retrieved in SetUp() still point to the
      // location of the sdMat data
      ASSERT_EQ(true,pointerUnchanged(mat1_,mat1_vp_,m_lda_));
      ASSERT_EQ(true,pointerUnchanged(mat2_,mat2_vp_,m_lda_));
      
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));

      // copy back in the changed block
      SUBR(sdMat_set_block)(mat1_,m1_copy,imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_SOUT(PHIST_DEBUG,"Changed block and copied it back in:\n");
      SUBR(sdMat_print)(mat1_,&iflag_);
      SUBR(sdMat_print)(m1_copy,&iflag_);
#endif      
      
      // check that the corresponding entries have changed
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(&mat1_vp_[MIDX(imin,jmin,m_lda_)],imax-imin+1,jmax-jmin+1,m_lda_,stride,val,mflag_));

      SUBR(sdMat_delete)(m1_copy,&iflag_);
      ASSERT_EQ(0,iflag_);
      }
    }

  TEST_F(CLASSNAME, sdMat_times_sdMat)
  {
    if (typeImplemented_ && nrows_ == ncols_ )
    {
      SUBR(sdMat_times_sdMat)(st::one(),mat1_,mat3_,st::one(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHIST_DEB("random*random+42");
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      SUBR(sdMat_from_device)(mat3_,&iflag_);

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols_; j++)
        {
          for(int k = 0; k < ncols_; k++)
          {
            mat2_vp_[MIDX(i,j,m_lda_)] -= mat1_vp_[MIDX(i,k,m_lda_)]*mat3_vp_[MIDX(k,j,m_lda_)];
          }
        }
      ASSERT_EQ(0,iflag_);
      // check result
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,(ST)42.0,mflag_),10*mt::eps());
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
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat3_,st::one(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHIST_DEB("random'*random+42");
#if PHIST_OUTLEV>=PHIST_DEBUG
/*
      SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
*/
#endif

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      SUBR(sdMat_from_device)(mat3_,&iflag_);

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
      {  
        for(int j = 0; j < ncols_; j++)
        {
          for(int k = 0; k < ncols_; k++)
          {
            mat2_vp_[MIDX(i,j,m_lda_)] -= 
            st::conj(mat1_vp_[MIDX(k,i,m_lda_)])*mat3_vp_[MIDX(k,j,m_lda_)];
          }
        }
      }
      // check result
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,(ST)42.0,mflag_),10*mt::eps());

    }
  }


  TEST_F(CLASSNAME, sdMat_times_sdMatT)
  {
    if (typeImplemented_ && nrows_ == ncols_ )
    {
#ifdef IS_COMPLEX
      // force some complex non-zero value even if complex_random fails
      mat1_vp_[0] = std::complex<MT>( st::imag(mat1_vp_[0]), st::real(mat1_vp_[0]) );
#endif
      SUBR(sdMat_times_sdMatT)(st::one(),mat1_,mat3_,st::one(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHIST_DEB("random'*random+42");
#if PHIST_OUTLEV>=PHIST_DEBUG
/*
      SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
*/
#endif

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      SUBR(sdMat_from_device)(mat3_,&iflag_);

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
      {  
        for(int j = 0; j < ncols_; j++)
        {
          for(int k = 0; k < ncols_; k++)
          {
            mat2_vp_[MIDX(i,j,m_lda_)] -= 
              mat1_vp_[MIDX(i,k,m_lda_)]*st::conj(mat3_vp_[MIDX(j,k,m_lda_)]);
          }
        }
      }
      // check result
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,(ST)42.0,mflag_),10*mt::eps());

    }
  }


  // check if we can construct an identity matrix
  TEST_F(CLASSNAME, identity)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      SUBR(sdMat_identity)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          ST val = st::zero();
          if( i == j )
            val = st::one();
          ASSERT_REAL_EQ(st::real(val),st::real(mat1_vp_[MIDX(i,j,m_lda_)]));
          ASSERT_REAL_EQ(st::imag(val),st::imag(mat1_vp_[MIDX(i,j,m_lda_)]));
        }
      }
    }
  }


  // the result of sdMat_random must be equal on all processes (if a comm. is given in sdMat_create!)
  TEST_F(CLASSNAME, parallel_random)
  {
    if( typeImplemented_ )
    {
      int stride = 1;

      int rank = 0;
      phist_comm_get_rank(comm_, &rank, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_random)(mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat1_vp_,nrows_,ncols_,m_lda_,stride,mflag_));

      // don't trick me (the test) by just using the same initilization for the random number generator on all processes!
      if( rank == 0 )
      {
        st::rand();
      }

      // do the same test again...
      SUBR(sdMat_random)(mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat1_vp_,nrows_,ncols_,m_lda_,stride,mflag_));
    }
  }

#if _NROWS_ >= 6 && _NCOLS_ >= 8
  TEST_F(CLASSNAME, happy_aliasing_views)
  {
    if( typeImplemented_ )
    {
      int stride = 1;
      SUBR(sdMat_random)(mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set up 4 block views in one matrix to multiply around
      TYPE(sdMat_ptr) v11 = NULL, v12 = NULL, v21 = NULL, v22 = NULL;
      SUBR(sdMat_view_block)(mat1_, &v11, 1,2, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &v21, 3,4, 1,3, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &v12, 0,2, 4,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &v22, 3,5, 5,6, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set up left/right and upper/lower block view to add around
      TYPE(sdMat_ptr) vu = NULL, vd = NULL, vl = NULL, vr = NULL;
      SUBR(sdMat_view_block)(mat1_, &vu, 0,2, 1,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &vd, 3,5, 1,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &vl, 0,4, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &vr, 1,5, 3,4, &iflag_);
      ASSERT_EQ(0,iflag_);


      // use "unaliased" other mat2_ as reference
      SUBR(sdMat_add_sdMat)(st::one(), mat1_, st::zero(), mat2_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set up 4 block views in one matrix to multiply around
      TYPE(sdMat_ptr) ref_v11 = NULL, ref_v12 = NULL, ref_v21 = NULL, ref_v22 = NULL;
      SUBR(sdMat_view_block)(mat2_, &ref_v11, 1,2, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_v21, 3,4, 1,3, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_v12, 0,2, 4,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_v22, 3,5, 5,6, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set up left/right and upper/lower block view to add around
      TYPE(sdMat_ptr) ref_vu = NULL, ref_vd = NULL, ref_vl = NULL, ref_vr = NULL;
      SUBR(sdMat_view_block)(mat2_, &ref_vu, 0,2, 1,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_vd, 3,5, 1,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_vl, 0,4, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_vr, 1,5, 3,4, &iflag_);
      ASSERT_EQ(0,iflag_);


      // do some operations, always start with result in reference mat2_ views than same on aliasing views in mat1_ and compare
      _ST_ alpha = 0.7;
      _ST_ beta = 0.3;

      // sdMat_add_sdMat up/down
      SUBR(sdMat_add_sdMat)(alpha, vu, beta, ref_vd, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(alpha, vu, beta, vd, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));

      // sdMat_add_sdMat left/right
      SUBR(sdMat_add_sdMat)(alpha, vr, beta, ref_vl, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(alpha, vr, beta, vl, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));

      // sdMat_times_sdMat 21,22->11
      SUBR(sdMat_times_sdMat)(alpha, v21, v22, beta, ref_v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_times_sdMat)(alpha, v21, v22, beta, v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));

      // sdMat_times_sdMat 22,21->12
      SUBR(sdMat_times_sdMat)(alpha, v22, v21, beta, ref_v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_times_sdMat)(alpha, v22, v21, beta, v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));


      // sdMatT_times_sdMat 22,22->11
      SUBR(sdMatT_times_sdMat)(alpha, v22, v22, beta, ref_v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMatT_times_sdMat)(alpha, v22, v22, beta, v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));

      // sdMatT_times_sdMat 21,21->12
      SUBR(sdMatT_times_sdMat)(alpha, v21, v21, beta, ref_v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMatT_times_sdMat)(alpha, v21, v21, beta, v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride));

      // delete ref views
      SUBR(sdMat_delete)(ref_vr, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_vl, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_vd, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_vu, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_v21, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(ref_v22, &iflag_);
      ASSERT_EQ(0,iflag_);

      // delete views
      SUBR(sdMat_delete)(vr, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(vl, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(vd, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(vu, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(v21, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(v22, &iflag_);
      ASSERT_EQ(0,iflag_);


    }
  }
#endif


  // check rank revealing cholesky decomposition
  TEST_F(CLASSNAME, cholesky)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      // -- check identity * 42 --
      SUBR(sdMat_identity)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(42.*st::one(), mat1_, st::zero(), mat2_, &iflag_);
      ASSERT_EQ(0,iflag_);
      // copy to mat1_
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      int rank = 0;
      int perm[nrows_];
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
//SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nrows_,rank);
      // assure that mat1_ is now permuted upper triangular
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i > j )
          {
            ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[MIDX(i,perm[j],m_lda_)]));
          }
        }
      }

      // check result
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));


      // -- check rank deficiency of last row/col --
      mat2_vp_[MIDX(nrows_-1,ncols_-1,m_lda_)] = st::zero();
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
//SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nrows_-1,rank);
      // assure that mat1_ is now permuted upper triangular
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i > j )
          {
            ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[MIDX(i,perm[j],m_lda_)]));
          }
        }
      }

      // check result
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));


      // -- create explicit hpd matrix from upper triangular part --
      SUBR(sdMat_put_value)(mat1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      int k = (nrows_*(nrows_+1))/2;
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i == j )
            mat1_vp_[MIDX(i,j,m_lda_)] = ST(10*k--);
          else if( i < j )
            mat1_vp_[MIDX(i,j,m_lda_)] = ST(k--);
          else
            mat1_vp_[MIDX(i,j,m_lda_)] = st::zero();
        }
      }
PHIST_SOUT(PHIST_INFO,"Predefined L^T:\n");
SUBR(sdMat_print)(mat1_,&iflag_);
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
PHIST_SOUT(PHIST_INFO,"M:\n");
SUBR(sdMat_print)(mat2_,&iflag_);
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
PHIST_SOUT(PHIST_INFO,"L^T:\n");
SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nrows_,rank);
      // assure that mat1_ is now permuted upper triangular
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i > j )
          {
            ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[MIDX(i,perm[j],m_lda_)]));
          }
        }
      }

      // check result
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(0,iflag_);
PHIST_SOUT(PHIST_INFO,"LL^T:\n");
SUBR(sdMat_print)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));


      // -- create explicit hermitian semi-positive definite matrix from upper triangular part --
      // requires working rank detection and pivoting!
      SUBR(sdMat_put_value)(mat1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      k = (nrows_*(nrows_-1))/2+nrows_-1;
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 1; j < ncols_; j++)
        {
          if( i <= j )
            mat1_vp_[MIDX(i,j,m_lda_)] = ST(k--);
          else
            mat1_vp_[MIDX(i,j,m_lda_)] = st::zero();
        }
      }
PHIST_SOUT(PHIST_INFO,"Predefined L^T:\n");
SUBR(sdMat_print)(mat1_,&iflag_);
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
PHIST_SOUT(PHIST_INFO,"M:\n");
SUBR(sdMat_print)(mat2_,&iflag_);
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
PHIST_SOUT(PHIST_INFO,"L^T:\n");
SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(nrows_-1,rank);
      // assure that mat1_ is now permuted upper triangular
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          if( i > j )
          {
            ASSERT_REAL_EQ(mt::zero(),st::abs(mat1_vp_[MIDX(i,perm[j],m_lda_)]));
          }
        }
      }

      // check result
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(0,iflag_);
PHIST_SOUT(PHIST_INFO,"LL^T:\n");
SUBR(sdMat_print)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat3_,mat2_),100*mt::eps());
    }
  }

