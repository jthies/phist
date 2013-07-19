#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithType< _ST_ >,
                 public KernelTestWithMap<_N_>
  {

public:

  typedef KernelTestWithVectors<_ST_,_N_,_M_> VTest;
  typedef KernelTestWithSdMats<_ST_,_M_,_M_> MTest;

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  
  //! V is n x m
  _TYPE_(mvec_ptr) V1_,V2_;

  //! M is m x m
  _TYPE_(sdMat_ptr) M1_,M2_;
  
  _ST_ *V1_vp_,*V2_vp_,*M1_vp_,*M2_vp_;
  
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  int ldaV1_,ldaV2_,ldaM1_,ldaM2_,stride_;
  
  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithType< _ST_ >::SetUp();
    KernelTestWithMap<_N_>::SetUp();
    if (this->typeImplemented_)
      {
      // create vectors V1 and V2, and vector views for setting/checking entries
      _SUBR_(mvec_create)(&V1_,this->map_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_extract_view)(V1_,&V1_vp_,&ldaV1_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_create)(&V2_,this->map_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(mvec_extract_view)(V2_,&V2_vp_,&ldaV2_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      // create matrices M1, M2 and views.
      _SUBR_(sdMat_create)(&M1_,this->m_,this->m_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(sdMat_extract_view)(M1_,&M1_vp_,&this->ldaM1_,&this->ierr_);
      _SUBR_(sdMat_create)(&M2_,this->m_,this->m_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      _SUBR_(sdMat_extract_view)(M2_,&M2_vp_,&this->ldaM2_,&this->ierr_);

      }
    stride_=1;
    }

  /*! Clean up.
   */
  virtual void TearDown()
    {
    if (this->typeImplemented_)
      {
      _SUBR_(mvec_delete)(V1_,&ierr_);
      _SUBR_(mvec_delete)(V2_,&ierr_);
      _SUBR_(sdMat_delete)(M1_,&ierr_);
      _SUBR_(sdMat_delete)(M2_,&ierr_);
      }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
    }

};

  // check if vectors are normalized correctly after QR factorization
  TEST_F(CLASSNAME, mvecT_times_mvec) 
    {
    if (typeImplemented_)
      {
      // fill V and W with ones
      _SUBR_(mvec_put_value)(V1_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      _SUBR_(mvec_put_value)(V2_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      _SUBR_(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      VTest::PrintVector(*cout,"ones",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      VTest::PrintVector(*cout,"ones",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"ones'*ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(M1_vp_,m_,m_,ldaM1_,stride_,(ST)nglob_));
      }
    }


