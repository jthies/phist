#ifndef CLASSNAME
#error "file not included correctly"
#endif
#if !defined(PHIST_HIGH_PRECISION_KERNELS) && defined(PHIST_HIGH_PRECISION_KERNELS_FORCE)
#define PHIST_HIGH_PRECISION_KERNELS
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_,_USE_VIEWS_> 
  {

public:

  typedef KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_,_USE_VIEWS_> MTest;

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
    MTest::TearDown();
  }

  /*! internal tests for forward/backward substition
   */
  void doForwardBackwardTestsWithPreparedMat3(int rank, int* perm);
};

  // check rank revealing cholesky decomposition
  TEST_F(CLASSNAME, cholesky)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      // -- check identity * 42 --
      SUBR(sdMat_identity)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_add_sdMat)((ST)42, mat1_, st::zero(), mat2_, &iflag_);
      ASSERT_EQ(0,iflag_);
      // copy to mat1_
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      int rank = 0;
      int perm[nrows_];
      int iflag_in=0;
#ifdef HIGH_PRECISION_KERNELS
      iflag_in=PHIST_ROBUST_REDUCTIONS;
#endif
      // these kernels only work on the host, so we need to manually perform up/downloads
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=iflag_in;
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
//SUBR(sdMat_print)(mat1_,&iflag_);
//      ASSERT_EQ(0,iflag_);
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
      iflag_=iflag_in;
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_HIGH_PRECISION_KERNELS
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));
#else
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat3_,mat2_),10*mt::eps());
#endif

      SUBR(sdMat_from_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // -- check rank deficiency of last row/col --
      mat2_vp_[MIDX(nrows_-1,ncols_-1,m_lda_)] = st::zero();

      SUBR(sdMat_to_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      rank = 0;
      iflag_=iflag_in;
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
//SUBR(sdMat_print)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nrows_-1,rank);
    SUBR(sdMat_to_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
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
      iflag_=iflag_in;
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_HIGH_PRECISION_KERNELS
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));
#else
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat3_,mat2_),10*mt::eps());
#endif

      // -- create explicit hpd matrix from upper triangular part --
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

      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      MTest::PrintSdMat(PHIST_DEBUG,"Predefined L^T",mat1_vp_,m_lda_,1,mpi_comm_);
      iflag_=iflag_in;
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
MTest::PrintSdMat(PHIST_DEBUG,"M",mat2_vp_,m_lda_,1,mpi_comm_);
      iflag_=iflag_in;
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      rank = 0;
      iflag_=iflag_in;
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
MTest::PrintSdMat(PHIST_DEBUG,"L^T",mat1_vp_,m_lda_,1,mpi_comm_);
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
      iflag_=iflag_in;
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      PrintSdMat(PHIST_DEBUG,"LL^T",mat3_vp_,m_lda_,1,mpi_comm_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_HIGH_PRECISION_KERNELS
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));
#else
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat3_,mat2_),10*mt::eps());
#endif


      // -- create explicit hermitian semi-positive definite matrix from upper triangular part --
      // requires working rank detection and pivoting!
      SUBR(sdMat_put_value)(mat1_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      k = (nrows_*(nrows_-1))/2+nrows_-1;
      for(int i = 0; i < nrows_; i++)
      {
        for(int j = 1; j < ncols_; j++)
        {
          if( i <= j )
            mat1_vp_[MIDX(i,j,m_lda_)] = ST(k--)+(ST)0.1*(ST)(j-i)*(ST)(mt::prand()-0.5)*st::cmplx_I();
          else
            mat1_vp_[MIDX(i,j,m_lda_)] = st::zero();
        }
      }

      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

PrintSdMat(PHIST_DEBUG,"predefined L^T",mat1_vp_,m_lda_,1,mpi_comm_);

      iflag_=iflag_in;
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"M",mat2_vp_,m_lda_,1,mpi_comm_);

      iflag_=iflag_in;
      SUBR(sdMat_add_sdMat)(st::one(), mat2_, st::zero(), mat1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // cholesky
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      rank = 0;
      iflag_=iflag_in;
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"L^T",mat1_vp_,m_lda_,1,mpi_comm_);

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
      iflag_=iflag_in;
      SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat1_,st::zero(),mat3_, &iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"LL^T",mat3_vp_,m_lda_,1,mpi_comm_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_HIGH_PRECISION_KERNELS
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(mat3_,mat2_));
#else
      ASSERT_NEAR(mt::one(),SdMatsEqual(mat3_,mat2_),10*mt::eps());
#endif
    }
  }

  // forward-backward substition
  TEST_F(CLASSNAME, forward_backward_subst)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      // data for substitution
      int rank = nrows_;
      int perm[nrows_];
      for(int i = 0; i < nrows_; i++)
        perm[i] = i;

      SUBR(sdMat_identity)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      PHIST_SOUT(PHIST_INFO, "Forward-backward substition test identity rhs\n");
      doForwardBackwardTestsWithPreparedMat3(rank,perm);
      if( HasFatalFailure() )
        return;

      // modify permutation
#if _NROWS_>2
        std::swap(perm[1],perm[0]);
#endif
#if _NROWS_ > 4
        std::swap(perm[3],perm[2]);
#endif
#if _NROWS_>5
        std::swap(perm[5],perm[1]);
#endif
      PHIST_SOUT(PHIST_INFO, "Forward-backward substition test with identity rhs and permutation\n");
      doForwardBackwardTestsWithPreparedMat3(rank,perm);
      if( HasFatalFailure() )
        return;

      for(int i = 0; i < nrows_; i++)
        perm[i] = i;

      SUBR(sdMat_random)(mat3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      PHIST_SOUT(PHIST_INFO, "Forward-backward substition test random rhs\n");
      doForwardBackwardTestsWithPreparedMat3(rank,perm);
      if( HasFatalFailure() )
        return;

      // modify permutation
#if (_NROWS_>2)
        std::swap(perm[1],perm[0]);
#endif
#if (_NROWS_>4)
        std::swap(perm[3],perm[2]);
#endif
#if (_NROWS_>5)
        std::swap(perm[5],perm[1]);
#endif
      PHIST_SOUT(PHIST_INFO, "Forward-backward substition test with random rhs and permutation\n");
      doForwardBackwardTestsWithPreparedMat3(rank,perm);
      if( HasFatalFailure() )
        return;
    }
  }

  void CLASSNAME::doForwardBackwardTestsWithPreparedMat3(int rank, int* perm)
  {
    int iflag_in=0;
#ifdef PHIST_HIGH_PRECISION_KERNELS
    iflag_in=PHIST_ROBUST_REDUCTIONS;
#endif
    // generate upper triangular factor ourselves
    SUBR(sdMat_put_value)(mat1_,st::zero(),&iflag_);
    ASSERT_EQ(0,iflag_);
    int k = (nrows_*(nrows_+1))/2;
    for(int i = 0; i < nrows_; i++)
    {
      for(int j = 0; j < ncols_; j++)
      {
        if( i <= j )
          mat1_vp_[MIDX(i,perm[j],m_lda_)] = ST(k--);
        else
          mat1_vp_[MIDX(i,perm[j],m_lda_)] = st::zero();
      }
    }
    
    SUBR(sdMat_to_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // multiply mat1_ with identity matrix with
PrintSdMat(PHIST_DEBUG,"X",mat3_vp_,m_lda_,1,mpi_comm_);
ASSERT_EQ(0,iflag_);
iflag_=iflag_in;
    SUBR(sdMat_times_sdMat)(st::one(),mat1_,mat3_,st::zero(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"R*X",mat2_vp_,m_lda_,1,mpi_comm_);

    // backward substitute
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    iflag_=iflag_in;
    SUBR(sdMat_backwardSubst_sdMat)(mat1_,perm,rank,mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_to_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"reconstructed X",mat2_vp_,m_lda_,1,mpi_comm_);
    // this should have reconstructed mat3_
iflag_=iflag_in;
    SUBR(sdMat_add_sdMat)(-st::one(),mat3_,st::one(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"difference",mat2_vp_,m_lda_,1,mpi_comm_);

ASSERT_EQ(0,iflag_);
#ifdef PHIST_HIGH_PRECISION_KERNELS
    ASSERT_NEAR(mt::one(),SdMatEqual(mat2_,st::zero()),100*mt::eps()*mt::eps());
#else
    ASSERT_NEAR(mt::one(),SdMatEqual(mat2_,st::zero()),10*mt::eps());
#endif
    // multiply mat1_^T with identity matrix with
    iflag_=iflag_in;
    SUBR(sdMatT_times_sdMat)(st::one(),mat1_,mat3_,st::zero(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"R^T*X",mat2_vp_,m_lda_,1,mpi_comm_);

    // forward substitute
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    iflag_=iflag_in;
    SUBR(sdMat_forwardSubst_sdMat)(mat1_,perm,rank,mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(sdMat_to_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(mvec_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
PrintSdMat(PHIST_DEBUG,"reconstructed X",mat2_vp_,m_lda_,1,mpi_comm_);

    // this should have reconstructed mat3_
    iflag_=iflag_in;
    SUBR(sdMat_add_sdMat)(-st::one(),mat3_,st::one(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
#ifdef PHIST_HIGH_PRECISION_KERNELS
    ASSERT_NEAR(mt::one(),SdMatEqual(mat2_,st::zero()),100*mt::eps()*mt::eps());
#else
    ASSERT_NEAR(mt::one(),SdMatEqual(mat2_,st::zero()),100*mt::eps());
#endif
  }

  // under construction... (TODO)
  TEST_F(CLASSNAME, DISABLED_pseudo_inverse)
  {
    if( typeImplemented_ )
    {
      int iflag_in=0;
#ifdef HIGH_PRECISION_KERNELS
      iflag_in=PHIST_ROBUST_REDUCTIONS;
#endif
      // copy mat2=mat1
      SUBR(sdMat_add_sdMat)(st::one(),mat1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // these kernels only work on the host, so we need to manually perform up/downloads
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // aliases for clarity: AplusT is A^{+,T}, the transpose of the Moore-Penrose Pseudo-Inverse
      TYPE(sdMat_ptr) A = mat2_, AplusT=mat1_;
      
      int rank;
      SUBR(sdMat_pseudo_inverse)(AplusT,&rank,&iflag_);
      ASSERT_EQ(0,iflag_);
      // for a random n x m matrix, we expect min(n,m) to be the rank
      ASSERT_EQ(std::min(nrows_,ncols_),rank);
      
      // check the four defining properties of A+:
#if 0      
      // 1. A A+ A = A
      TripleProduct(A,false,AplusT,true,A,false,mat3_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(st::one(),SdMatsEqual(A,mat3_),st::eps());
      
      // 2. A+ A A+ = A+ => A+^T A^T A+^T = A+^T
      TripleProduct(AplusT,false,A,false,AplusT,false,mat3_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(st::one(),SdMatsEqual(AplusT,mat3_),st::eps());
#endif
      // 3. (AA+)^* = AA+
      TYPE(sdMat_ptr) mat_tmp=NULL;
      SUBR(sdMat_create)(&mat_tmp, nrows_,nrows_,comm_,&iflag_);
      SUBR(sdMat_times_sdMatT)(st::one(),A,AplusT,st::zero(),mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(mt::one(),mt::one()+MTest::symmetry_check(mat_tmp,&iflag_),mt::eps());
      SUBR(sdMat_delete)(mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // 4. (A+ A)^* = A+ A
      SUBR(sdMat_create)(&mat_tmp, ncols_,ncols_,comm_,&iflag_);
      SUBR(sdMatT_times_sdMat)(st::one(),AplusT,A,st::zero(),mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(mt::one(),mt::one()+MTest::symmetry_check(mat_tmp,&iflag_),mt::eps());
      SUBR(sdMat_delete)(mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
