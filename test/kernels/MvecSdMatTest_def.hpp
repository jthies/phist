#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithType< _ST_ >,
                 public KernelTestWithMap<_N_>
{

public:

  typedef KernelTestWithVectors<_ST_,_N_,_M_> V1Test;
  typedef KernelTestWithVectors<_ST_,_N_,_K_> V2Test;
  typedef KernelTestWithSdMats<_ST_,_M_,_M_> MTest;

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  //TODO: k is currently ignored in this test
  static const int k_=_K_;
  
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
      SUBR(mvec_create)(&V1_,this->map_,this->m_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_put_value)(V1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_extract_view)(V1_,&V1_vp_,&ldaV1_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_create)(&V2_,this->map_,this->k_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(mvec_put_value)(V2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_extract_view)(V2_,&V2_vp_,&ldaV2_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      // create matrices M1, M2 and views.
      SUBR(sdMat_create)(&M1_,this->m_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(M1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(M1_,&M1_vp_,&this->ldaM1_,&this->iflag_);
      SUBR(sdMat_create)(&M2_,this->m_,this->k_,this->comm_,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_put_value)(M2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_extract_view)(M2_,&M2_vp_,&this->ldaM2_,&this->iflag_);
    }
    stride_=1;
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    if (this->typeImplemented_)
    {
      SUBR(mvec_delete)(V1_,&iflag_);
      SUBR(mvec_delete)(V2_,&iflag_);
      SUBR(sdMat_delete)(M1_,&iflag_);
      SUBR(sdMat_delete)(M2_,&iflag_);
    }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
  }

  /*! compare sdMats on several procs
   */
  void SUBR(sdMat_parallel_check_)(TYPE(const_sdMat_ptr) mat, int* iflag)
  {
    *iflag = 0;
    // TODO: use correct communicator
    int n,m;
    PHIST_CHK_IERR(SUBR(sdMat_from_device)((void*)mat,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &m, iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &n, iflag),*iflag);
    _ST_* buff = new _ST_[m*n];
    _ST_* mat_raw;
    lidx_t lda;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)((TYPE(sdMat_ptr))mat, &mat_raw, &lda, iflag),*iflag);
    // copy data to buffer
    for(int j = 0; j < n; j++)
      for(int i = 0; i < m; i++)
        buff[j*m+i] = mat_raw[MIDX(i,j,lda)];
    // broadcast
    PHIST_CHK_IERR(*iflag = MPI_Bcast(buff,m*n,::phist::ScalarTraits<_ST_>::mpi_type(),0,MPI_COMM_WORLD),*iflag);
    // check
    int error = 0;
    for(int j = 0; j < n; j++)
      for(int i = 0; i < m; i++)
        if( buff[j*m+i] != mat_raw[MIDX(i,j,lda)] )
          error = 1;
    int globError = 0;
    PHIST_CHK_IERR(*iflag = MPI_Allreduce(&error,&globError,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD),*iflag);

    delete[] buff;

    if( globError )
    {
#if PHIST_OUTLEV>=PHIST_DEBUG
      PHIST_CHK_IERR(SUBR(sdMat_print)(mat,iflag),*iflag);
#endif
      *iflag = -1;
      return;
    }
  }
};

  // check ones(n,m)'*ones(n,k)=n*ones(m,k)
  TEST_F(CLASSNAME, mvecT_times_mvec) 
  {
    if (typeImplemented_)
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(V2_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      V1Test::PrintVector(*cout,"ones",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      V2Test::PrintVector(*cout,"ones",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"ones'*ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
      ASSERT_REAL_EQ(mt::one(),SdMatEqual(M1_,(ST)nglob_));
      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check ones(n,m)'*ones(n,m)=n*ones(m,m)
  TEST_F(CLASSNAME, mvecT_times_mvec_with_manual_comparison) 
  {
    if (typeImplemented_)
    {
      // fill V and W with random numbers
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
    
      MT err_r=mt::zero();
      MT err_i=mt::zero();
      for (int i=0; i<m_; i++)
        for (int j=0; j<k_; j++)
        {
          ST err_ij=st::zero();
          // manually compute V1(:,i)*V2(:,j) using compensated summation ("Kahan")
          // for achieving good accuracy
          // (http://en.wikipedia.org/wiki/Kahan_summation_algorithm)
          ST c = st::zero(); // a running compensation for lost low-order bits
          for (int ii=0; ii< nloc_; ii++)
          {
            ST y = st::conj(V1_vp_[VIDX(ii,i,ldaV1_)])
                 *          V2_vp_[VIDX(ii,j,ldaV2_)] - c; // c is 0 in the beginning
            ST t = err_ij + y;          // err_ij is big, y small, so low-order digits of y are lost.
            c = (t - err_ij) - y;  // (t - err_ij) recovers the high-order part of y; 
                                   // subtracting y recovers -(low part of y)
            err_ij = t;           // Algebraically, c should always 
                                  // be zero. Beware overly-aggressive optimizing 
                                  // compilers!
                                  // Next time around, the lost low part will 
                                  // be added to y in a fresh attempt.

          }
#ifdef HAVE_MPI
          // TODO - we should use Kahan summation here as well, but since our tests are not
          // carried out on many processes, it doesn't matter right now.
          ST err_ij_local=err_ij;
          MPI_Allreduce(&err_ij_local,&err_ij,1,st::mpi_type(),MPI_SUM,mpi_comm_);
#endif
          err_ij-=M1_vp_[MIDX(i,j,ldaM1_)];
          err_r=std::max(err_r,st::real(err_ij));
          err_i=std::max(err_i,st::imag(err_ij));
        }
      ASSERT_NEAR(mt::zero(),err_r,100*mt::eps());
      ASSERT_NEAR(mt::zero(),err_i,100*mt::eps());
    }
  }

  // check ones(n,m)*ones(m,k)=m*ones(n,k)
  TEST_F(CLASSNAME, mvec_times_sdMat)
  {
    if (typeImplemented_)
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(V2_,(MT)42.0*st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      V1Test::PrintVector(*cout,"ones",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
      V2Test::PrintVector(*cout,"ones*ones",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
#endif
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V2_,(ST)m_));
      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

#if (_M_==_K_)

  // check ones(n,m)*ones(m,m)=m*ones(n,m)
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place)
  {
    if (typeImplemented_)
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat_inplace)(V1_,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      MTest::PrintSdMat(*cout,"ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
      V1Test::PrintVector(*cout,"ones*ones",V1_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
#endif
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V1_,(ST)m_));
      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check ones(n,m)*ones(m,m)=m*ones(n,m)
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place_with_random_data)
  {
    if (typeImplemented_)
    {
      // fill V and W with ones
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat_inplace)(V1_,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      V1Test::PrintVector(*cout,"result_inplace",V1_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
      V2Test::PrintVector(*cout,"result_out_of_place",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
#endif
      ASSERT_REAL_EQ(mt::one(),MvecsEqual(V1_,V2_));
      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check that we can zero out some columns of V by multiplying
  // a view of them with a zero sdMat.
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place_with_views)
  {
    if (typeImplemented_)
    {
      // fill V and M with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // need to view a square block from M1
      int bs=std::min(std::min(k_,m_),3);
      
      // create a view of a part of V and M
      int imin=0;
      int jmin=0;
      
      if (m_>bs)
      {
        imin++;
      }
      
      if (m_>bs)
      {
        jmin++;
      }

      int imax=std::min(m_-1, imin+bs-1);
      int jmax=std::min(m_-1, jmin+bs-1);
      
      ASSERT_EQ(jmax-jmin,imax-imin);

      mvec_ptr_t V=NULL;
      SUBR(mvec_view_block)(V1_,&V, jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      sdMat_ptr_t M=NULL;
      SUBR(sdMat_view_block)(M1_,&M, imin,imax,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_put_value)(M,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHIST_SOUT(PHIST_DEBUG,"range of zero M-block: (%d:%d,%d:%d)",imin,imax,jmin,jmax);
      V1Test::PrintVector(*cout,"1-vec",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"1-mat with hole",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif

      ASSERT_REAL_EQ(mt::one(),MvecEqual(V,st::one()));
      
      SUBR(mvec_times_sdMat_inplace)(V,M,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),MvecEqual(V,st::zero()));

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHIST_SOUT(PHIST_DEBUG,"range of zero M-block: (%d:%d,%d:%d)",imin,imax,jmin,jmax);
      V1Test::PrintVector(*cout,"1-vec with hole",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"1-mat with hole",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
            
      SUBR(mvec_from_device)(V, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check all vector entries, something like [1 .. 1 0 .. 0 1 .. 1] in all rows
      ASSERT_REAL_EQ(mt::one(),ArrayEqual( &V1_vp_[VIDX(0,0,ldaV1_)],
        nloc_,jmin,ldaV1_,stride_,st::one(),vflag_));
      ASSERT_REAL_EQ(mt::one(),ArrayEqual( &V1_vp_[VIDX(0,jmin,ldaV1_)],
        nloc_,jmax-jmin+1,ldaV1_,stride_,st::zero(),vflag_));
      ASSERT_REAL_EQ(mt::one(),ArrayEqual( &V1_vp_[VIDX(0,jmax+1,ldaV1_)],
        nloc_,m_-jmax-1,ldaV1_,stride_,st::one(),vflag_));

    }
  }
  /* m==k */
#endif

  // random check
  TEST_F(CLASSNAME, random_mvecT_times_mvec) 
  {
    if (typeImplemented_)
    {
      // fill V and W with ones
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      V1Test::PrintVector(*cout,"random",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      V2Test::PrintVector(*cout,"random",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"random'*random",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
      SUBR(sdMat_parallel_check_)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
      // to check the result, scale V1 by -1, add V1'V2 again and compare with 0
      SUBR(mvec_scale)(V1_,-st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::one(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
  
      // TODO - how much accuracy should we demand? With Kahan summation it should
      //        be accurate to machine precision, but without it probably won't...
      ASSERT_NEAR(mt::one(),SdMatEqual(M1_,st::zero()),100*mt::eps());
    
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
        if( off1[i]+m1[i] > m_ || off2[i]+m2[i] > k_  || off1_M[i]+m1[i] > m_ || off2_M[i]+m2[i] > k_)
          continue;
        PHIST_SOUT(PHIST_DEBUG, "Test offsets: off1: %d, off2: %d, m1: %d, m2: %d, off1_M: %d, off2_M: %d\n", off1[i], off2[i], m1[i], m2[i], off1_M[i], off2_M[i]);

        // create views to parts of mvecs and sdmats
        mvec_ptr_t V1 = NULL;
        SUBR(mvec_view_block)(V1_,&V1,off1[i],off1[i]+m1[i]-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        mvec_ptr_t V2 = NULL;
        SUBR(mvec_view_block)(V2_,&V2,off2[i],off2[i]+m2[i]-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        sdMat_ptr_t M1 = NULL;
        SUBR(sdMat_view_block)(M1_,&M1,off1_M[i],off1_M[i]+m1[i]-1,off2_M[i],off2_M[i]+m2[i]-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        sdMat_ptr_t M2 = NULL;
        SUBR(sdMat_view_block)(M2_,&M2,off1[i],off1[i]+m1[i]-1,off2[i],off2[i]+m2[i]-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        sdMat_ptr_t M3=NULL;
        SUBR(sdMat_create)(&M3, m1[i], m2[i], comm_, &iflag_);
        ASSERT_EQ(0,iflag_);
        _ST_* M3_vp;
        lidx_t lda_M3;
        SUBR(sdMat_extract_view)(M3,&M3_vp,&lda_M3, &iflag_);
        ASSERT_EQ(0,iflag_);

        PHIST_DEB("Note: we are just using views inside the random vectors\n");
        PHIST_DEB("col-range V1: [%d:%d]\n",off1[i],off1[i]+m1[i]-1);
        PHIST_DEB("col-range V2: [%d:%d]\n",off2[i],off2[i]+m2[i]-1);
        PHIST_DEB("idx-range M:  [%d:%d,%d:%d]\n",off1_M[i],off1_M[i]+m1[i]-1,off2_M[i],off2_M[i]+m2[i]-1);

        // set V1 and V2 to 0,
        // fill (viewed) V and W with random numbers
        SUBR(mvec_put_value)(V1_,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_random)(V1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_put_value)(V2_,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_random)(V2,&iflag_);
        ASSERT_EQ(0,iflag_);
        // fill M1_ and M2_ with zeros
        SUBR(sdMat_put_value)(M1_,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_put_value)(M2_,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_put_value)(M3,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);

        // first compute the full matrix V'W, giving non-zeros
        // only in a block with offset (off1, off2) and size m1 x m2
        SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M2_,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(sdMat_from_device)(M2_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // viewed vectors into regular sdMat
        SUBR(mvecT_times_mvec)(st::one(),V1,V2,st::zero(),M3,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(sdMat_from_device)(M3,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        // now the version with views
        SUBR(mvecT_times_mvec)(st::one(),V1,V2,st::zero(),M1,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(sdMat_from_device)(M1_,&iflag_);
        ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
        SUBR(mvec_from_device)(V1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_from_device)(V2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        V1Test::PrintVector(*cout,"random",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
        V2Test::PrintVector(*cout,"random",V2_vp_,nloc_,ldaV2_,stride_,mpi_comm_);
        
                MTest::PrintSdMat(*cout,"random'*random without views",M2_vp_,ldaM2_,stride_,mpi_comm_);
        PHIST_SOUT(PHIST_DEBUG,"viewed block in result");
        SUBR(sdMat_print)(M2,&iflag_);
        ASSERT_EQ(0,iflag_);

        PHIST_SOUT(PHIST_DEBUG,"view'*view, non-viewed result");
        SUBR(sdMat_print)(M3,&iflag_);
        ASSERT_EQ(0,iflag_);

        std::ostringstream ss;
        ss<<"rnd'*rnd in location ("<<off1_M[i]<<":"<<off1_M[i]+m1[i]-1<<","
                                         <<off2_M[i]<<":"<<off2_M[i]+m2[i]-1<<")";
        MTest::PrintSdMat(*cout,ss.str().c_str(), M1_vp_,ldaM1_,stride_,mpi_comm_);
        PHIST_SOUT(PHIST_DEBUG,"viewed block in result");
        SUBR(sdMat_print)(M1,&iflag_);
        ASSERT_EQ(0,iflag_);
#endif
        SUBR(sdMat_parallel_check_)(M1_,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(sdMat_parallel_check_)(M2_,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(sdMat_parallel_check_)(M3,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        // subtract the result without views from the one with non-viewed target sdMat
        SUBR(sdMat_add_sdMat)(-st::one(),M1, st::one(),M3,&iflag_);
        ASSERT_EQ(0,iflag_);

        // the result should be zero!
        SUBR(sdMat_from_device)(M3,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),ArrayEqual(M3_vp,m1[i],m2[i],lda_M3,stride_,st::zero(),mflag_),200*mt::eps());

        // subtract the two viewed blocks in the result sdMats
        SUBR(sdMat_add_sdMat)(-st::one(),M1, st::one(),M2,&iflag_);
        ASSERT_EQ(0,iflag_);

        // the result should be zero!
        SUBR(sdMat_from_device)(M2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),ArrayEqual(M2_vp_,m_,k_,ldaM2_,stride_,st::zero(),mflag_),200*mt::eps());

        // clean up at the end of the loop
        SUBR(mvec_delete)(V1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(V2,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_delete)(M1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_delete)(M2,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sdMat_delete)(M3,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }
  }

const int CLASSNAME::k_;
const int CLASSNAME::m_;
const int CLASSNAME::n_;
