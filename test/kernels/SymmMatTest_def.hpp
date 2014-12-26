#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_> 
{

  public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    
    if (typeImplemented_)
      {
      SUBR(read_mat)("sprandsym",nglob_,&A1_,&iflag_);
      if ( A1_==NULL )
        {
        haveMats_=false;
        }
      else
        {
        haveMats_=true;
        }
      }
    }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    if (typeImplemented_)
    {
      ASSERT_EQ(0,delete_mat(A1_));
    }
  }

// the matrices may have individual maps, so we need to recreate all vectors with the specific map of the matrix!
void rebuildVectors(TYPE(const_sparseMat_ptr) A)
{
  if (typeImplemented_ && haveMats_)
  {
    // set vec1 to be a valid X, vec2 and vec3 a valid Y in Y=AX
    const_map_ptr_t range_map, domain_map;
    SUBR(sparseMat_get_range_map)(A,&range_map,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sparseMat_get_domain_map)(A,&domain_map,&iflag_);

    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);

    lidx_t lda;
    SUBR(mvec_create)(&vec1_,domain_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec1_,&vec1_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    SUBR(mvec_create)(&vec2_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec2_,&vec2_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    SUBR(mvec_create)(&vec3_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec3_,&vec3_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    phist_map_get_local_length(domain_map, &nloc_, &iflag_);
    ASSERT_EQ(0,iflag_);
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
    if (typeImplemented_)
    {
      ASSERT_TRUE(AssertNotNull(A1_));
    }
  }

  TEST_F(CLASSNAME, A1_probe_symmetry)
  {
    if (typeImplemented_ && haveMats_)
    {

      // matrices may have different maps
      rebuildVectors(A1_);

      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A1_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
      TYPE(sdMat_ptr) VtAV=NULL;
      SUBR(sdMat_create)(&VtAV,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),vec1_,vec2_,st::zero(),VtAV,&iflag_);
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
      ASSERT_NEAR(mt::zero(),diff,1000*mt::eps());
      ASSERT_TRUE(sum>1000*mt::eps());

      SUBR(sdMat_delete)(VtAV,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
