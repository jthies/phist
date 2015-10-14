#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_> 
{

  public:
    typedef KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_> VTest;

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    VTest::SetUp();

    haveMats_ = false;
    if (typeImplemented_ && !problemTooSmall_)
    {
      if( nglob_ == 20 )
      {
        SUBR(read_mat)("symmMat",comm_,nglob_,&A1_,&iflag_);
      }
      else
      {
        SUBR(read_mat)("sprandsym",comm_,nglob_,&A1_,&iflag_);
      }
      ASSERT_EQ(0,iflag_);
      if( A1_ != NULL )
      {
        haveMats_ = true;

        const_map_ptr_t range_map;
        SUBR(sparseMat_get_range_map)(A1_,&range_map,&iflag_);
        replaceMap(range_map);
      }
    }
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    VTest::TearDown();
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_EQ(0,delete_mat(A1_));
    }
  }


TYPE(sparseMat_ptr) A1_; // some random symmetric matrix

protected:

int delete_mat(TYPE(sparseMat_ptr) A)
{
  if (A!=NULL)
  {
    SUBR(sparseMat_delete)(A,&iflag_);
  }
  return iflag_;
}

  bool haveMats_;
};

  TEST_F(CLASSNAME, read_matrices) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_TRUE(AssertNotNull(A1_));
    }
  }

  TEST_F(CLASSNAME, A1_probe_symmetry)
  {
    if (typeImplemented_ && !problemTooSmall_ && haveMats_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A1_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
      TYPE(sdMat_ptr) VtAV=NULL;
      SUBR(sdMat_create)(&VtAV,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),vec1_,vec2_,st::zero(),VtAV,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(VtAV,&iflag_);
#endif
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(VtAV,&iflag_);
      ASSERT_EQ(0,iflag_);

      _ST_* C_raw=NULL;
      int lda;
      SUBR(sdMat_extract_view)(VtAV,&C_raw,&lda,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check V'AV != 0 and V'AV symmetric      
      _MT_ sum=mt::one();;
      _MT_ diff=mt::zero();
      for (int i=0;i<_NV_;i++)
        for (int j=i;j<_NV_;j++)
        {
          _ST_ cij=C_raw[i*lda+j];
          _ST_ cji=C_raw[j*lda+i];
          sum=std::min(sum,st::abs(cij)+st::abs(cji));
          diff=std::max(diff,st::abs(st::conj(cij)-cji));
        }
      ASSERT_NEAR(mt::zero(),diff,sqrt(mt::eps()));
      ASSERT_TRUE(sum>sqrt(mt::eps()));

      SUBR(sdMat_delete)(VtAV,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
