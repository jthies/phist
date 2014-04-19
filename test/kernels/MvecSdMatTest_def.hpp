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
  TYPE(mvec_ptr) V1_,V2_;

  //! M is m x m
  TYPE(sdMat_ptr) M1_,M2_;
  
  _ST_ *V1_vp_,*V2_vp_,*M1_vp_,*M2_vp_;
  
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  lidx_t ldaV1_,ldaV2_,ldaM1_,ldaM2_,stride_;
  
  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithType< _ST_ >::SetUp();
    KernelTestWithMap<_N_>::SetUp();
    if (this->typeImplemented_)
      {
      // create vectors V1 and V2, and vector views for setting/checking entries
      SUBR(mvec_create)(&V1_,this->map_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_put_value)(V1_,st::zero(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_extract_view)(V1_,&V1_vp_,&ldaV1_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_create)(&V2_,this->map_,this->m_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(mvec_put_value)(V2_,st::zero(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_extract_view)(V2_,&V2_vp_,&ldaV2_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      // create matrices M1, M2 and views.
      SUBR(sdMat_create)(&M1_,this->m_,this->m_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(sdMat_put_value)(M1_,st::zero(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(sdMat_extract_view)(M1_,&M1_vp_,&this->ldaM1_,&this->ierr_);
      SUBR(sdMat_create)(&M2_,this->m_,this->m_,this->comm_,&this->ierr_);
      ASSERT_EQ(0,this->ierr_);
      SUBR(sdMat_put_value)(M2_,st::zero(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(sdMat_extract_view)(M2_,&M2_vp_,&this->ldaM2_,&this->ierr_);
      }
    stride_=1;
    }

  /*! Clean up.
   */
  virtual void TearDown()
    {
    if (this->typeImplemented_)
      {
      SUBR(mvec_delete)(V1_,&ierr_);
      SUBR(mvec_delete)(V2_,&ierr_);
      SUBR(sdMat_delete)(M1_,&ierr_);
      SUBR(sdMat_delete)(M2_,&ierr_);
      }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
    }

  /*! compare sdMats on several procs
   */
  void SUBR(sdMat_parallel_check_)(TYPE(const_sdMat_ptr) mat, int* ierr)
  {
    *ierr = 0;
    // TODO: use correct communicator
    int n,m;
    PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &m, ierr),*ierr);
    PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &n, ierr),*ierr);
    _ST_* buff = new _ST_[m*n];
    _ST_* mat_raw;
    lidx_t lda;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)((TYPE(sdMat_ptr))mat, &mat_raw, &lda, ierr),*ierr);
    // copy data to buffer
    for(int j = 0; j < n; j++)
      for(int i = 0; i < m; i++)
        buff[j*m+i] = mat_raw[MIDX(i,j,lda)];
    // broadcast
    PHIST_CHK_IERR(*ierr = MPI_Bcast(buff,m*n,::phist::ScalarTraits<_ST_>::mpi_type(),0,MPI_COMM_WORLD),*ierr);
    // check
    int error = 0;
    for(int j = 0; j < n; j++)
      for(int i = 0; i < m; i++)
        if( buff[j*m+i] != mat_raw[MIDX(i,j,lda)] )
          error = 1;
    int globError = 0;
    PHIST_CHK_IERR(*ierr = MPI_Allreduce(&error,&globError,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD),*ierr);

    delete[] buff;

    if( globError )
    {
      PHIST_CHK_IERR(SUBR(sdMat_print)(mat,ierr),*ierr);
      *ierr = -1;
      return;
    }
  }
};

  // check ones(n,m)'*ones(n,m)=n*ones(m,m)
  TEST_F(CLASSNAME, mvecT_times_mvec) 
    {
    if (typeImplemented_)
      {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_put_value)(V2_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      VTest::PrintVector(*cout,"ones",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      VTest::PrintVector(*cout,"ones",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"ones'*ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(M1_vp_,m_,m_,ldaM1_,stride_,(ST)nglob_,mflag_));
      SUBR(sdMat_parallel_check_)(M1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      }
    }

  // check ones(n,m)*ones(m,m)=m*ones(n,m)
  TEST_F(CLASSNAME, mvec_times_sdMat) 
    {
    if (typeImplemented_)
      {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_put_value)(V2_,(MT)42.0*st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(sdMat_put_value)(M1_,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,st::zero(),V2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      VTest::PrintVector(*cout,"ones",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
      VTest::PrintVector(*cout,"ones*ones",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
#endif
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(V2_vp_,nloc_,m_,ldaV2_,stride_,(ST)m_,vflag_));
      SUBR(sdMat_parallel_check_)(M1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      }
    }

  // random check
  TEST_F(CLASSNAME, random_mvecT_times_mvec) 
    {
    if (typeImplemented_)
      {
      // fill V and W with ones
      SUBR(mvec_random)(V1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_random)(V2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      VTest::PrintVector(*cout,"random",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      VTest::PrintVector(*cout,"random",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"random'*random",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
      SUBR(sdMat_parallel_check_)(M1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      }
    }

  // random check with partial views of partial mvecs and sdMats
  TEST_F(CLASSNAME, random_mvecT_times_mvec_with_inside_views) 
    {
    if (typeImplemented_ && m_ > 4)
      {
      std::vector<int> off1;
      std::vector<int> off2;
      std::vector<int> m1;
      std::vector<int> m2;
      std::vector<int> off1_M;
      std::vector<int> off2_M;

      off1.push_back(0);  m1.push_back(3);  off2.push_back(1);  m2.push_back(3);  off1_M.push_back(0);  off2_M.push_back(1);
      off1.push_back(1);  m1.push_back(1);  off2.push_back(0);  m2.push_back(6);  off1_M.push_back(4);  off2_M.push_back(0);
      off1.push_back(3);  m1.push_back(5);  off2.push_back(0);  m2.push_back(1);  off1_M.push_back(2);  off2_M.push_back(3);
      off1.push_back(7);  m1.push_back(2);  off2.push_back(4);  m2.push_back(3);  off1_M.push_back(1);  off2_M.push_back(2);

      for(int i = 0; i < off1.size(); i++)
      {
        if( off1[i]+m1[i] > m_ || off2[i]+m2[i] > m_  || off1_M[i]+m1[i] > m_ || off2_M[i]+m2[i] > m_)
          continue;

        // create views to parts of mvecs and sdmats
        mvec_ptr_t V1 = NULL;
        SUBR(mvec_view_block)(V1_,&V1,off1[i],off1[i]+m1[i]-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        mvec_ptr_t V2 = NULL;
        SUBR(mvec_view_block)(V2_,&V2,off2[i],off2[i]+m2[i]-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        sdMat_ptr_t M1 = NULL;
        SUBR(sdMat_view_block)(M1_,&M1,off1_M[i],off1_M[i]+m1[i]-1,off2_M[i],off2_M[i]+m2[i]-1,&ierr_);
        ASSERT_EQ(0,ierr_);

        // fill V and W with ones
        SUBR(mvec_put_value)(V1_,st::zero(),&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_random)(V1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_put_value)(V2_,st::zero(),&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_random)(V2,&ierr_);
        ASSERT_EQ(0,ierr_);
        // fill M1_ with zeros
        SUBR(sdMat_put_value)(M1_,st::zero(),&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvecT_times_mvec)(-st::one(),V1,V2,st::one(),M1,&ierr_);
        ASSERT_EQ(0,ierr_);

        PHIST_DEB("Note: we are just using views inside the random vectors");
#if PHIST_OUTLEV>=PHIST_DEBUG
        VTest::PrintVector(*cout,"random",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
        VTest::PrintVector(*cout,"random",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
        MTest::PrintSdMat(*cout,"zero-random'*random",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
        SUBR(sdMat_parallel_check_)(M1_,&ierr_);
        ASSERT_EQ(0,ierr_);

        PHIST_DEB("check the result with full mvecT_times_mvec");
        // copy M1 back to correct position
        SUBR(sdMat_put_value)(M2_,st::zero(),&ierr_);
        ASSERT_EQ(0,ierr_);
        sdMat_ptr_t M2 = NULL;
        SUBR(sdMat_view_block)(M2_,&M2,off1[i],off1[i]+m1[i]-1,off2[i],off2[i]+m2[i]-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(sdMat_add_sdMat)(st::one(),M1,st::zero(),M2,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(sdMat_delete)(M2,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::one(),M2_,&ierr_);
        ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
        MTest::PrintSdMat(*cout,"zero-random'*random+random'*random",M2_vp_,ldaM2_,stride_,mpi_comm_);
#endif
        SUBR(sdMat_parallel_check_)(M2_,&ierr_);
        ASSERT_EQ(0,ierr_);

        // the result should be zero!
        ASSERT_NEAR(mt::one(),ArrayEqual(M2_vp_,m_,m_,ldaM2_,stride_,st::zero()),200*mt::eps());

        SUBR(mvec_delete)(V1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(V2,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(sdMat_delete)(M1,&ierr_);
        ASSERT_EQ(0,ierr_);
      }
      }
    }
