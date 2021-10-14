/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef CLASSNAME
#error "file not included correctly"
#endif
#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif

#ifdef PHIST_KERNEL_LIB_TPETRA__disabled
#define DO_TPETRA_TESTS 1
#endif

#ifdef DO_TPETRA_TESTS
#include "Tpetra_MultiVector.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_View.hpp"
#include "Teuchos_FancyOStream.hpp"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_,_USE_VIEWS_> 
{

public:

  typedef KernelTestWithSdMats<ST,_NROWS_,_NCOLS_,_USE_VIEWS_> MTest;

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    MTest::SetUp();
    if( typeImplemented_ )
    {
      SUBR(sdMat_random)(mat1_,&iflag_);
        ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(mat2_,ST(42.0),&iflag_);
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
    MTest::TearDown();
  }

#if DO_TPETRA_TESTS

  bool check_entries(TYPE(sdMat_ptr) M, int imin, int imax, int jmin, int jmax)
  {
/*
    std::cout << "range: "<<imin<<":"<<imax<<","<<jmin<<":"<<jmax<<std::endl;
    Teuchos::FancyOStream ofs(Teuchos::rcp(&std::cout,false));
    ((Tpetra::MultiVector<ST,phist_lidx,phist_gidx>*)M)->describe(ofs,Teuchos::VERB_EXTREME);
*/
    ST* M_raw=nullptr;
    phist_lidx ldM;
    int nrows, ncols;
    int iflag;
    SUBR(sdMat_extract_view)(M,&M_raw, &ldM, &iflag);
    if (iflag) return false;
    SUBR(sdMat_get_nrows)(M,&nrows,&iflag);
    if (iflag) throw iflag;
    SUBR(sdMat_get_ncols)(M,&ncols,&iflag);
    if (iflag) throw iflag;
    if (nrows!=(imax-imin+1)) return false;
    if (ncols!=(jmax-jmin+1)) return false;
    bool result=true;
    for (int j=0; j<ncols; j++)
      for (int i=0; i<nrows; i++)
      {
        int ii=imin+i+1;
        int jj=jmin+j+1;
        result &= M_raw[MIDX(i,j,ldM)]==ST(ii*1000+jj);
      }
#if PHIST_OUTLEV>=PHIST_DEBUG
    if (!result) 
    {
      std::ostringstream os;
      os << "Tpetra view test: M("<<imin<<":"<<imax<<","<<jmin<<":"<<jmax<<") looks incorrect:\n";
      for (int i=0; i<nrows; i++)
      {
        for (int j=0; j<ncols; j++)
        {
          os << "\t"<<M_raw[MIDX(i,j,ldM)];
        }
        os << "\n";
      }
      std::cout << os.str()<<std::endl;
    }
#endif
    return result;
  }
  
  bool check_entries_recursive(TYPE(sdMat_ptr) M, int imin, int imax, int jmin, int jmax)
  {
    TYPE(sdMat_ptr) Msub=nullptr;
    int iflag=0;
    bool result_this_level=check_entries(M,imin,imax,jmin,jmax);
    if (!result_this_level || (imax-imin==1 && jmax-jmin==1)) return result_this_level;
    bool other_results=true;
    for (int i0=imin; i0<=imax; i0++)
    for (int i1=i0; i1<=imax; i1++)
    for (int j0=jmin; j0<=jmax; j0++)
    for (int j1=j0; j1<=jmax; j1++)
    {
      if (i0==imin && i1==imax && j0==jmin && j1==jmax) continue;
      int nrsub=i1-i0+1, ncsub=j1-j0+1;
      SUBR(sdMat_view_block)(M,&Msub,i0-imin,i1-imin,j0-jmin,j1-jmin,&iflag);
      if (iflag) throw iflag;
      bool this_result=check_entries_recursive(Msub,i0,i1,j0,j1);
      if (this_result==false) 
      {
        other_results=false;
        PHIST_OUT(PHIST_DEBUG,"subview(%d:%d,%d:%d) was incorrect!\n",i0,i1,j0,j1);
        break;
      }
    }
    if (Msub!=nullptr) SUBR(sdMat_delete)(Msub,&iflag);
    if (iflag) throw iflag;
    return other_results;
  }
#endif
  /*! internal tests for forward/backward substition
   */
  void doForwardBackwardTestsWithPreparedMat3(int rank, int* perm);
};


#ifdef DO_TPETRA_TESTS

  TEST_F(CLASSNAME, tpetra_data_layout)
  {
    if (!typeImplemented_) return;
    typedef Tpetra::MultiVector<ST,phist_lidx,phist_gidx> tpetra_mvec;
    ASSERT_TRUE(m_lda_>=nrows_);
    tpetra_mvec* M=(tpetra_mvec*)mat1_;
    auto Mview = M->getLocalView<Kokkos::HostSpace>();
    bool offsets_match_lda=true;
    for (int j=0; j<ncols_-1; j++)
    {
      for (int i=0; i<nrows_; i++)
      {
        offsets_match_lda &= (&(Mview(i,j+1)) - &(Mview(i,j))==m_lda_);
      }
    }
    ASSERT_TRUE(offsets_match_lda);
  }

  TEST_F(CLASSNAME, DISABLED_tpetra_nested_view_block)
  {
    typedef Tpetra::MultiVector<ST,phist_lidx,phist_gidx> tpetra_mvec;

    tpetra_mvec* M=(tpetra_mvec*)mat1_;
    auto Mkokkos = M->getLocalView<Kokkos::HostSpace>();
    
    for (int i=0; i<nrows_; i++)
    {
      for (int j=0; j<ncols_; j++)
      {
        Mkokkos(i,j) = ST((i+1)*1000+(j+1));
      }
    }
    bool all_views_correct=check_entries_recursive(mat1_,0,nrows_-1,0,ncols_-1);
    ASSERT_TRUE(all_views_correct);
  }

#endif


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
    // and do not crash. Note that it is not easy to verify the  
    // correctness of the functions in our current settings: the 
    // semantics are different if there is an accelerator, and we
    // do not have a function to test wether ther is one.
    if (typeImplemented_)
    {
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  TEST_F(CLASSNAME, put_value_done_on_host_and_device)
  {
    if (typeImplemented_)
    {
      int stride=1;
      ST val=st::prand();
      SUBR(sdMat_put_value)(mat1_,val,&iflag_);
      ASSERT_EQ(0,iflag_);
      EXPECT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_,nrows_,ncols_,m_lda_,stride,val,mflag_));
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      EXPECT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_,nrows_,ncols_,m_lda_,stride,val,mflag_));
    }
  }

  TEST_F(CLASSNAME, random_same_on_host_and_device)
  {
    if (typeImplemented_)
    {
      SUBR(sdMat_random)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ST sum1=st::zero();
      for (int i=0; i<nrows_; i++)
      {
        for (int j=0; j<ncols_; j++)
        {
          sum1 += mat1_vp_[MIDX(i,j,m_lda_)];
        }
      }
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ST sum2=st::zero();
      for (int i=0; i<nrows_; i++)
      {
        for (int j=0; j<ncols_; j++)
        {
          sum2 += mat1_vp_[MIDX(i,j,m_lda_)];
        }
      }
      ASSERT_REAL_EQ(st::real(sum1),st::real(sum2));
      ASSERT_REAL_EQ(st::imag(sum1),st::imag(sum2));
      ASSERT_TRUE(sum1!=st::zero());
    }
  }

  TEST_F(CLASSNAME, identity_done_on_host_and_device)
  {
    if (typeImplemented_)
    {
      int stride=1;
      SUBR(sdMat_identity)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_identity)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int i=0; i<std::min(nrows_,ncols_); i++)
      {
        mat1_vp_[i*m_lda_+i] -= st::one();
        mat2_vp_[i*m_lda_+i] -= st::one();
      }
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(mat1_vp_,nrows_,ncols_,m_lda_,stride,st::zero(),mflag_));
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,stride,st::zero(),mflag_));
    }
  }

  TEST_F(CLASSNAME, sync_values)
  {
    if (typeImplemented_)
    {
      // manually set values depending on MPI rank
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          mat1_vp_[MIDX(i,j,m_lda_)] = (mpi_rank_*nrows_+i)*ncols_+j;
          mat2_vp_[MIDX(i,j,m_lda_)] = i*ncols_+j; // this should be the result on all procs
        }
      }
      
      // synchronize
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_to_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(mat1_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // make sure everyone has the values set by rank 0
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
    }
  }

  TEST_F(CLASSNAME, sync_values_of_view)
  {
    if (typeImplemented_)
    {
      ST a0=ST((mpi_rank_+1)*1234);
      SUBR(sdMat_put_value)(mat1_,a0,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(mat2_,a0,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      int i0=std::min(1,nrows_-1);
      int i1=nrows_-1;
      int j0=std::min(1,ncols_-1);
      int j1=ncols_-1;
      // manually set values depending on MPI rank
      for(int i = i0; i <= i1; i++)
      {
        for(int j = j0; j <= j1; j++)
        {
          mat1_vp_[MIDX(i,j,m_lda_)] = (mpi_rank_*nrows_+i)*ncols_+j;
          mat2_vp_[MIDX(i,j,m_lda_)] = i*ncols_+j; // this should be the result on all procs
        }
      }
      
      // synchronize
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_to_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      TYPE(sdMat_ptr) mat1v=nullptr;
      SUBR(sdMat_view_block)(mat1_,&mat1v,i0,i1,j0,j1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(mat1v,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1v,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(mat1v,&iflag_);
      ASSERT_EQ(0,iflag_);
      // make sure everyone has the values set by rank 0
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
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
      ASSERT_EQ(true,pointerUnchanged(mat1_,mat1_vp_,m_lda_));
      ASSERT_EQ(true,pointerUnchanged(mat2_,mat2_vp_,m_lda_));
    }

  }

  // X = 0*Y + a*X = a*X
  TEST_F(CLASSNAME, scale_by_axpy)
  {
    if (typeImplemented_)
    {
      ST alpha = st::zero();
      ST beta  = st::prand();
      PHIST_OUT(9,"axpy, alpha=%f+%f i, beta=%f+%f i\n",st::real(alpha),
          st::imag(alpha),st::real(beta),st::imag(beta));
#if PHIST_OUTLEV>=PHIST_DEBUG
      MTest::PrintSdMat(PHIST_DEBUG,"before scale",mat2_vp_,m_lda_,1,mpi_comm_);
#endif
      SUBR(sdMat_add_sdMat)(alpha,mat1_,beta,mat2_,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      MTest::PrintSdMat(PHIST_DEBUG,"after scale",mat2_vp_,m_lda_,1,mpi_comm_);
#endif
      ASSERT_NEAR(mt::one(),SdMatEqual(mat2_,beta*ST(42.0)),mt::sqrt(mt::eps()));
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
          mat1_vp_[MIDX(i,j,m_lda_)] = alpha*mat1_vp_[MIDX(i,j,m_lda_)]+ST(42.0)*beta;
        }
      }
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat1_,mat2_),mt::sqrt(mt::eps()));
    }
  }

  // X = a*Y' + b*X
  TEST_F(CLASSNAME, random_sdMatT_add_sdMat)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      ST beta = st::prand();
      ST alpha= st::prand();

      SUBR(sdMatT_add_sdMat)(alpha,mat1_,beta,mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // calculate solution by hand
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 0; j < ncols_; j++)
        {
          mat3_vp_[MIDX(i,j,m_lda_)] = alpha*st::conj(mat1_vp_[MIDX(j,i,m_lda_)])+ST(42.0)*beta;
        }
      }
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat3_,mat2_),mt::sqrt(mt::eps()));
    }
  }


  // view certain rows and columns of a serial dense matrix,
  // manipulate them and check the original matrix has changed.
  TEST_F(CLASSNAME, view_block)
  {
    if (typeImplemented_)
    {
      int imin=std::min(2,nrows_-1);
      int imax=std::min(4,nrows_-1);
      int jmin=std::min(1,ncols_-1);
      int jmax=std::min(3,ncols_-1);
      int stride=1;
      TYPE(sdMat_ptr) m1_view=nullptr;
      SUBR(sdMat_view_block)(mat1_,&m1_view,imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      ST* val_ptr;
      phist_lidx lda;
      SUBR(sdMat_extract_view)(m1_view,&val_ptr,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      if( jmax-jmin > 1 ) 
      {
        ASSERT_EQ(lda,m_lda_);
      }

      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+MIDX(imin,jmin,m_lda_),val_ptr,imax-imin+1,jmax-jmin+1,m_lda_,stride));

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(mat1_,&iflag_);
      SUBR(sdMat_print)(m1_view,&iflag_);
#endif
      
      SUBR(sdMat_extract_view)(m1_view,&val_ptr,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      if( jmax-jmin > 1 ) 
      {
        ASSERT_EQ(lda,m_lda_);
      }
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+MIDX(imin,jmin,m_lda_),val_ptr,imax-imin+1,jmax-jmin+1,m_lda_,stride,mflag_));

      // set all the viewed entries to a certain value and check that the original vector is
      // changed.
      ST val = random_number();
      SUBR(sdMat_put_value)(m1_view,val,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),ArraysEqual(mat1_vp_+MIDX(imin,jmin,m_lda_),val_ptr,imax-imin+1,jmax-jmin+1,m_lda_,stride));
      
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_print)(m1_view,&iflag_);
#endif
      
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(&mat1_vp_[MIDX(imin,jmin,m_lda_)],imax-imin+1,jmax-jmin+1,m_lda_,stride,val,mflag_));

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
      ST outer_val = st::prand();
      SUBR(sdMat_put_value)(mat1_,outer_val,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // now create a view
      int imin=std::min(3,nrows_-1);
      int imax=std::min(7,nrows_-1);
      int jmin=std::min(2,ncols_-1);
      int jmax=std::min(5,ncols_-1);
      TYPE(sdMat_ptr) view = nullptr;
      SUBR(sdMat_view_block)(mat1_,&view,imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      // set the data in the view to some other value
      ST view_val = st::prand();
      SUBR(sdMat_put_value)(view,view_val,&iflag_);
      ASSERT_EQ(0,iflag_);

      // view part of view
      int nrows_view = imax-imin+1;
      int ncols_view = jmax-jmin+1;
      ASSERT_EQ(0,iflag_);
      int imin2=std::min(2,nrows_view-1);
      int imax2=std::min(4,nrows_view-1);
      int jmin2=std::min(1,ncols_view-1);
      int jmax2=std::min(3,ncols_view-1);
      TYPE(sdMat_ptr) view2 = nullptr;
      SUBR(sdMat_view_block)(view, &view2, imin2,imax2,jmin2, jmax2, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      // set data in the inner view to yet another value
      ST inner_val = st::prand();
      SUBR(sdMat_put_value)(view2, inner_val, &iflag_);
      ASSERT_EQ(0,iflag_);

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
      TYPE(sdMat_ptr) m1_copy=nullptr;
      
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
      ST* val_ptr;
      phist_lidx lda;
      SUBR(sdMat_extract_view)(m1_copy,&val_ptr,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      //note: can't use ArraysEqual here because we (may) have different strides      
      MT maxval=mt::zero();
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
      ST val = random_number();
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
      EXPECT_EQ(0,iflag_);
      SUBR(sdMat_print)(m1_copy,&iflag_);
      EXPECT_EQ(0,iflag_);
#endif      
      
      // check that the corresponding entries have changed
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(&mat1_vp_[MIDX(imin,jmin,m_lda_)],imax-imin+1,jmax-jmin+1,m_lda_,stride,val,mflag_));

      SUBR(sdMat_delete)(m1_copy,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  TEST_F(CLASSNAME, get_set_data_CM)
  {
    if (typeImplemented_)
    {
      phist_lidx data_lda = nrows_+ncols_+3;
      _ST_ data_rm[nrows_*data_lda], data_cm[ncols_*data_lda];

      // test get_data from col-major to col-major. mat1_ is random right now
      SUBR(sdMat_get_data)(mat1_,data_cm,data_lda,0,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(mt::one(),ArraysEqualWithDifferentLDA(mat1_vp_,data_cm,nrows_,ncols_,m_lda_,data_lda,1,0,mt::eps()));

      SUBR(sdMat_set_data)(mat2_,data_cm,data_lda,0,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
    }
  }

  TEST_F(CLASSNAME, get_set_data_RM)
  {
    if (typeImplemented_)
    {
      phist_lidx data_lda = nrows_+ncols_+3;
      _ST_ data_rm[nrows_*data_lda], data_cm[ncols_*data_lda];
      // test get_data from col-major to row-major. mat1_ is random right now
      SUBR(sdMat_get_data)(mat1_,data_rm,data_lda,1,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int i=0; i<nrows_; i++)
      {
        for (int j=0; j<ncols_; j++)
        {
          data_cm[j*data_lda+i]=data_rm[i*data_lda+j];
        }
      }
      ASSERT_EQ(mt::one(),ArraysEqualWithDifferentLDA(mat1_vp_,data_cm,nrows_,ncols_,m_lda_,data_lda,1,0,mt::eps()));

      SUBR(sdMat_set_data)(mat2_,data_rm,data_lda,1,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat1_,mat2_));
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
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,ST(42.0),mflag_),mt::sqrt(mt::eps()));
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
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,ST(42.0),mflag_),mt::sqrt(mt::eps()));

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
      ASSERT_NEAR(mt::one(),ArrayEqual(mat2_vp_,nrows_,ncols_,m_lda_,1,ST(42.0),mflag_),mt::sqrt(mt::eps()));

    }
  }

  TEST_F(CLASSNAME, products_involving_nonsquare_sdMats)
  {
    if (typeImplemented_ )
    {
      phist_lidx ncols2=std::max(3,ncols_-1);
      TYPE(sdMat_ptr) matA=mat1_, matB=nullptr, matC=nullptr;
      SUBR(sdMat_create)(&matB,ncols_,ncols2,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_create)(&matC,nrows_,ncols2,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      phist_lidx ldA=m_lda_,ldB,ldC;
      ST *matA_vp=mat1_vp_,*matB_vp=nullptr,*matC_vp=nullptr;
      
      SUBR(sdMat_extract_view)(matB,&matB_vp,&ldB,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(matC,&matC_vp,&ldC,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // A is (nrows  x ncols)
      // B is (ncols  x ncols2)
      // C is (nrows  x ncols2)

      //////////////////////////////////////////
      // 1) C = C + A*B                       //
      //////////////////////////////////////////

      SUBR(sdMat_random)(matA,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(matB,&iflag_);
      ASSERT_EQ(0,iflag_);
      ST valC=st::one()*ST(42.0) - st::cmplx_I()*ST(3.14);
      SUBR(sdMat_put_value)(matC,valC,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(sdMat_times_sdMat)(st::one(),matA,matB,st::one(),matC,&iflag_);
      ASSERT_EQ(0,iflag_);

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols2; j++)
        {
          for(int k = 0; k < ncols_; k++)
          {
            matC_vp[MIDX(i,j,ldC)] -= matA_vp[MIDX(i,k,ldA)]*matB_vp[MIDX(k,j,ldB)];
          }
        }
      // check result
      ASSERT_NEAR(mt::one(),SdMatEqual(matC,valC),mt::sqrt(mt::eps()));

      //////////////////////////////////////////
      // 2) B = B + A^H*C                     //
      //////////////////////////////////////////

      SUBR(sdMat_random)(matA,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(matC,&iflag_);
      ASSERT_EQ(0,iflag_);
      ST valB=st::one()*ST(1.234) - st::cmplx_I()*(ST(2.99));
      SUBR(sdMat_put_value)(matB,valB,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMatT_times_sdMat)(st::one(),matA,matC,st::one(),matB,&iflag_);
      ASSERT_EQ(0,iflag_);

      // subtract matrix product by hand
      for(int i = 0; i < ncols_; i++)
        for(int j = 0; j < ncols2; j++)
        {
          for(int k = 0; k < nrows_; k++)
          {
            matB_vp[MIDX(i,j,ldB)] -= st::conj(matA_vp[MIDX(k,i,ldA)])*matC_vp[MIDX(k,j,ldC)];
          }
        }
      // check result
      ASSERT_NEAR(mt::one(),SdMatEqual(matB,valB),mt::sqrt(mt::eps()));

      //////////////////////////////////////////
      // 2) A = A + C*B^H                     //
      //////////////////////////////////////////

      SUBR(sdMat_random)(matC,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(matB,&iflag_);
      ASSERT_EQ(0,iflag_);
      ST valA=st::one()*ST(1.5) - st::cmplx_I()*ST(5.1);
      SUBR(sdMat_put_value)(matA,valA,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(sdMat_times_sdMatT)(st::one(),matC,matB,st::one(),matA,&iflag_);
      ASSERT_EQ(0,iflag_);

      // subtract matrix product by hand
      for(int i = 0; i < nrows_; i++)
        for(int j = 0; j < ncols_; j++)
        {
          for(int k = 0; k < ncols2; k++)
          {
            matA_vp[MIDX(i,j,ldA)] -= matC_vp[MIDX(i,k,ldC)]*st::conj(matB_vp[MIDX(j,k,ldB)]);
          }
        }
      // check result
      ASSERT_NEAR(mt::one(),SdMatEqual(matA,valA),mt::sqrt(mt::eps()));

      // clean up
      SUBR(sdMat_delete)(matB,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(matC,&iflag_);
      ASSERT_EQ(0,iflag_);

    }
  }

  // check if we can construct an identity matrix
  TEST_F(CLASSNAME, identity)
  {
    if( typeImplemented_ )
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

#ifndef PHIST_KERNEL_LIB_BUILTIN
      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat1_vp_,nrows_,ncols_,m_lda_,stride,mflag_));

      // don't trick me (the test) by just using the same initilization for the random number generator on all processes!
      if( rank == 0 )
      {
        st::rand();
      }

      // do the same test again...
      SUBR(sdMat_random)(mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      
#ifndef PHIST_KERNEL_LIB_BUILTIN
      SUBR(sdMat_sync_values)(mat1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat1_vp_,nrows_,ncols_,m_lda_,stride,mflag_));

#ifdef PHIST_KERNEL_LIB_BUILTIN
      // verify prand doesn't "outsync" sdMat_random
      ST val = st::prand();
      SUBR(sdMat_random)(mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat1_vp_,nrows_,ncols_,m_lda_,stride,mflag_));
#endif
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
      TYPE(sdMat_ptr) v11 = nullptr, v12 = nullptr, v21 = nullptr, v22 = nullptr;
      SUBR(sdMat_view_block)(mat1_, &v11, 1,2, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &v21, 3,4, 1,3, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &v12, 0,2, 4,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat1_, &v22, 3,5, 5,6, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set up left/right and upper/lower block view to add around
      TYPE(sdMat_ptr) vu = nullptr, vd = nullptr, vl = nullptr, vr = nullptr;
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
      TYPE(sdMat_ptr) ref_v11 = nullptr, ref_v12 = nullptr, ref_v21 = nullptr, ref_v22 = nullptr;
      SUBR(sdMat_view_block)(mat2_, &ref_v11, 1,2, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_v21, 3,4, 1,3, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_v12, 0,2, 4,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_v22, 3,5, 5,6, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set up left/right and upper/lower block view to add around
      TYPE(sdMat_ptr) ref_vu = nullptr, ref_vd = nullptr, ref_vl = nullptr, ref_vr = nullptr;
      SUBR(sdMat_view_block)(mat2_, &ref_vu, 0,2, 1,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_vd, 3,5, 1,6, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_vl, 0,4, 1,2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_view_block)(mat2_, &ref_vr, 1,5, 3,4, &iflag_);
      ASSERT_EQ(0,iflag_);


      // do some operations, always start with result in reference mat2_ views than same on aliasing views in mat1_ and compare
      ST alpha = 0.7;
      ST beta = 0.3;

      // sdMat_add_sdMat up/down
      SUBR(sdMat_add_sdMat)(alpha, vu, beta, ref_vd, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(alpha, vu, beta, vd, &iflag_);
      ASSERT_EQ(0,iflag_);
     
      ASSERT_NEAR   (mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride),mt::sqrt(mt::eps()));

      // sdMat_add_sdMat left/right
      SUBR(sdMat_add_sdMat)(alpha, vr, beta, ref_vl, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)(alpha, vr, beta, vl, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR   (mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride),mt::sqrt(mt::eps()));

      // sdMat_times_sdMat 21,22->11
      SUBR(sdMat_times_sdMat)(alpha, v21, v22, beta, ref_v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_times_sdMat)(alpha, v21, v22, beta, v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR   (mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride),mt::sqrt(mt::eps()));

      // sdMat_times_sdMat 22,21->12
      SUBR(sdMat_times_sdMat)(alpha, v22, v21, beta, ref_v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_times_sdMat)(alpha, v22, v21, beta, v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR   (mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride),mt::sqrt(mt::eps()));


      // sdMatT_times_sdMat 22,22->11
      SUBR(sdMatT_times_sdMat)(alpha, v22, v22, beta, ref_v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMatT_times_sdMat)(alpha, v22, v22, beta, v11, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR   (mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride),mt::sqrt(mt::eps()));

      // sdMatT_times_sdMat 21,21->12
      SUBR(sdMatT_times_sdMat)(alpha, v21, v21, beta, ref_v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMatT_times_sdMat)(alpha, v21, v21, beta, v12, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR   (mt::one(),ArraysEqual(mat1_vp_,mat2_vp_,nrows_,ncols_,m_lda_,stride),mt::sqrt(mt::eps()));

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

