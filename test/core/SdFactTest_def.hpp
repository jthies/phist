/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
  
  //! given mat2_, with expected rank, run some tests for the qb transofrmation
  //! (core of the SVQB algorithm by Stathopoulos)
  void run_qb_test(int expected_rank)
  {
      // create A as a product B'B to make sure it is spd
      SUBR(sdMatT_times_sdMat)(st::one(),mat2_,mat2_,st::zero(),mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(sdMat_from_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      MT nrms_ref[nrows_];
      for (int i=0; i<nrows_; i++) nrms_ref[i]=std::sqrt(st::abs(mat1_vp_[i*m_lda_+i]));

      // qb transform: A->B^ such that Q=V*B^ would normalize the original V
      int rank = -42;
      int iflag_in=0;
#ifdef HIGH_PRECISION_KERNELS
      iflag_in=PHIST_ROBUST_REDUCTIONS;
#endif

SUBR(sdMat_print)(mat1_,&iflag_);
ASSERT_EQ(0,iflag_);

      iflag_=iflag_in;
      _MT_ rankTol=mt::rankTol(iflag_in==PHIST_ROBUST_REDUCTIONS);
      _MT_ nrmsV[nrows_];
      SUBR(sdMat_qb)(mat1_,mat3_,nrmsV,&rank,rankTol,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_to_device)(mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);

SUBR(sdMat_print)(mat1_,&iflag_);
ASSERT_EQ(0,iflag_);
      ASSERT_EQ(expected_rank,rank);
      // check that the squareroots of the diagonal elements are correctly returned
      for(int i = 0; i < nrows_; i++)
      {
        ASSERT_REAL_EQ(nrmsV[i],nrms_ref[i]);
      }

      // check that the inverse is correctly returned. In the case of a rank-deficient input
      // matrix we get that Bi*B is a 'rank identity' matrix [I 0; 0 0] with the last n-r rows
      // and columns zero.
      SUBR(sdMat_put_value)(mat2_,ST(0),&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int i=0; i<rank; i++) mat2_vp_[i*m_lda_+i]=ST(1);
      SUBR(sdMat_to_device)(mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(mat2_,&iflag_);
      iflag_=iflag_in;
      SUBR(sdMat_times_sdMat)(-st::one(),mat3_,mat1_,st::one(),mat2_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_print)(mat1_,&iflag_);
      SUBR(sdMat_print)(mat3_,&iflag_);
      SUBR(sdMat_print)(mat2_,&iflag_);
      // note that we will get some zeros on the diagonal if the matrix doesn't have full rank,
      // hence the non-standard test for equality below.
#ifdef PHIST_HIGH_PRECISION_KERNELS
      ASSERT_REAL_EQ(MT(1),SdMatEqual(mat2_,st::zero()));
#else
      ASSERT_NEAR(MT(1),SdMatEqual(mat2_,st::zero()),10*mt::eps());
#endif
  }


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

#if (_N_==_M_)
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
      _MT_ rankTol=mt::rankTol(iflag_in==PHIST_ROBUST_REDUCTIONS);
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,rankTol,&iflag_);
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
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,rankTol,&iflag_);
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
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,rankTol,&iflag_);
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
            mat1_vp_[MIDX(i,j,m_lda_)] = ST(k--)+ST(0.1*(j-i))*ST(mt::prand()-0.5)*st::cmplx_I();
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
      SUBR(sdMat_cholesky)(mat1_,perm,&rank,rankTol,&iflag_);
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

  TEST_F(CLASSNAME, qb_random_matrix)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      // -- check with full rank=m
      SUBR(sdMat_random)(mat2_,&iflag_);
      run_qb_test(nrows_);
    }
  }

  // -- check with rank=1
  TEST_F(CLASSNAME, qb_constant_matrix)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      SUBR(sdMat_put_value)(mat2_,ST(1),&iflag_);
      run_qb_test(1);
    }
  }

  // -- check with rank=1
  TEST_F(CLASSNAME, qb_zero_matrix)
  {
    if( typeImplemented_ && nrows_ == ncols_ )
    {
      SUBR(sdMat_put_value)(mat2_,ST(0),&iflag_);
      run_qb_test(0);
    }
  }

#endif /* _N_==_M_ */

  TEST_F(CLASSNAME, svd)
  {
    if( !typeImplemented_ ) return;
    // for random A, compute A=U*Sigma*V' and verify the product gives A,
    // Sigma is diagonal and it's entries are decreasing and positive, and
    // U and V are orthonormal.
    TYPE(sdMat_ptr) A=mat1_, A_bak=mat2_, Sigma=mat3_, U=NULL, Vt=NULL,UtU=NULL,VtV=NULL;

    SUBR(sdMat_create)(&U,nrows_,nrows_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_create)(&Vt,ncols_,ncols_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_create)(&UtU,nrows_,nrows_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_create)(&VtV,ncols_,ncols_,comm_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
#ifdef IS_DOUBLE
    _MT_ tol=mt::eps()*100;
#else
    _MT_ tol=mt::eps()*100;
#endif

    SdMatOwner<_ST_> _U(U),_Vt(Vt),_UtU(UtU),_VtV(VtV);

    // initialize Sigma, U and V with random entries, this should not hurt...
    //SUBR(sdMat_put_value)(Sigma,st::zero(),&iflag_);
    SUBR(sdMat_random)(Sigma,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_random)(U,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_random)(Vt,&iflag_);
    ASSERT_EQ(0,iflag_);
            
    //save A
    SUBR(sdMat_add_sdMat)(st::one(),A,st::zero(),A_bak,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    // compute the singular value decomposition, destroying A
    SUBR(sdMat_svd)(A,U,Sigma,Vt,&iflag_);
    ASSERT_EQ(0,iflag_);

    // see if we can reconstruct A:
    triple_product(U,false,Sigma,false,Vt,false,A,&iflag_);
    ASSERT_EQ(0,iflag_);

    ASSERT_NEAR(mt::one(),SdMatsEqual(A,A_bak),nrows_*ncols_*tol);
    _ST_* Sigma_vp=NULL;
    phist_lidx ldSigma;
    SUBR(sdMat_extract_view)(Sigma,&Sigma_vp,&ldSigma,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check singular values are positive and decreasing on the diagonal
    _ST_ last_sigma=st::zero();
    for (int i=std::min(nrows_,ncols_)-1; i>=0; i--)
    {
      _ST_ sigma_i=Sigma_vp[MIDX(i,i,ldSigma)];
      ASSERT_TRUE(st::real(sigma_i)>=st::real(last_sigma));
      ASSERT_TRUE(mt::abs(st::imag(sigma_i))<=mt::eps());
      last_sigma=sigma_i;
    }
    _MT_ max_offdiag=mt::zero();
    for (int i=0; i<nrows_; i++)
    {
      for (int j=0; j<ncols_; j++)
      {
        if (i!=j) max_offdiag=std::max(max_offdiag,st::abs(Sigma_vp[MIDX(i,j,ldSigma)]));
      }
    }
    ASSERT_NEAR(mt::zero(),max_offdiag,mt::eps());
    // compute U'U and V'V and check they are identity matrices
    SUBR(sdMat_identity)(UtU,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMatT_times_sdMat)(-st::one(),U,U,st::one(),UtU,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),SdMatEqual(UtU,mt::zero()),10*mt::eps());

    SUBR(sdMat_identity)(VtV,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_times_sdMatT)(-st::one(),Vt,Vt,st::one(),VtV,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),SdMatEqual(VtV,mt::zero()),100*mt::eps());
  }

  TEST_F(CLASSNAME, pseudo_inverse)
  {
    if( typeImplemented_ )
    {
#ifdef IS_DOUBLE
      _MT_ tol=mt::eps()*100;
#else
      // we get quite poor accuracy in SP, probably the better thing would be to
      // store sdMats in DP in that case!
      _MT_ tol=0.00002;
#endif
      int iflag_in=0;
#ifdef HIGH_PRECISION_KERNELS
      iflag_in=PHIST_ROBUST_REDUCTIONS;
#endif
      // copy mat2=mat1
      SUBR(sdMat_add_sdMat)(st::one(),mat1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // aliases for clarity: AplusT is A^{+,T}, the transpose of the Moore-Penrose Pseudo-Inverse
      TYPE(sdMat_ptr) A = mat2_, AplusT=mat1_, A_chk=mat3_;

      int rank;
      _MT_ rankTol=mt::rankTol(iflag_in==PHIST_ROBUST_REDUCTIONS);
      iflag_=iflag_in;
      SUBR(sdMat_pseudo_inverse)(AplusT,&rank,rankTol,&iflag_);
      ASSERT_EQ(0,iflag_);
      // for a random n x m matrix, we expect min(n,m) to be the rank
      ASSERT_EQ(std::min(nrows_,ncols_),rank);

      // check the four defining properties of A+:

      // 1. A A+ A = A
      MTest::triple_product(A,false,AplusT,true,A,false,A_chk,&iflag_);
      ASSERT_EQ(0,iflag_);
      EXPECT_NEAR(mt::one(),SdMatsEqual(A,A_chk),nrows_*nrows_*tol);
      
      // 2. A+ A A+ = A+ => A+^T A^T A+^T = A+^T
      MTest::triple_product(AplusT,false,A,true,AplusT,false,A_chk,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      EXPECT_NEAR(mt::one(),SdMatsEqual(AplusT,A_chk),ncols_*ncols_*tol);

      // 3. (AA+)^* = AA+
      TYPE(sdMat_ptr) mat_tmp=NULL;
      SUBR(sdMat_create)(&mat_tmp, nrows_,nrows_,comm_,&iflag_);
      SUBR(sdMat_times_sdMatT)(st::one(),A,AplusT,st::zero(),mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);

      EXPECT_NEAR(mt::one(),mt::one()+MTest::symmetry_check(mat_tmp,&iflag_),5*nrows_*nrows_*mt::eps());
      SUBR(sdMat_delete)(mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // 4. (A+ A)^H = A+ A
      mat_tmp=NULL;
      SUBR(sdMat_create)(&mat_tmp, ncols_,ncols_,comm_,&iflag_);
      SUBR(sdMatT_times_sdMat)(st::one(),AplusT,A,st::zero(),mat_tmp,&iflag_);
      EXPECT_EQ(0,iflag_);

      EXPECT_NEAR(mt::one(),mt::one()+MTest::symmetry_check(mat_tmp,&iflag_),100*mt::eps());
      SUBR(sdMat_delete)(mat_tmp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
