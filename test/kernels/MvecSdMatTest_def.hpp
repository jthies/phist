#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly"
#endif

#ifdef USE_VIEWS
#undef USE_VIEWS
#endif
#if _USE_VIEWS_V1_ || _USE_VIEWS_V2_ || _USE_VIEWS_M_
#define USE_VIEWS 1
#else
#define USE_VIEWS 0
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithType< _ST_ >,
                 public virtual KernelTestWithMap<_N_>
{

public:

  class V1Test : public KernelTestWithVectors<_ST_,_N_,_M_,_USE_VIEWS_V1_> {
    public: virtual void TestBody(){}
  };
  class V2Test : public KernelTestWithVectors<_ST_,_N_,_K_,_USE_VIEWS_V2_> {
    public: virtual void TestBody(){}
  };
  class MTest : public KernelTestWithSdMats<_ST_,_M_,_K_,_USE_VIEWS_M_> {
    public: virtual void TestBody(){}
  };

  //! mvec/sdMat sizes
  static const int n_=_N_;
  static const int m_=_M_;
  //TODO: k is currently ignored in this test
  static const int k_=_K_;
  
  //! V is n x m
  TYPE(mvec_ptr) V1_ = NULL, V2_ = NULL;

  //! M is m x m
  TYPE(sdMat_ptr) M1_ = NULL, M2_ = NULL;
  
  _ST_ *V1_vp_,*V2_vp_,*M1_vp_,*M2_vp_;
  
  // how defines the data layout. Vector
  // i starts at (i-1)*lda. Entries j and j+1
  // are at memory locations (i-1)*lda+stride*j
  // and (i-1)*lda+stride*(j+1), respectively.
  lidx_t ldaV1_,ldaV2_,ldaM1_,ldaM2_,stride_;

  V1Test v1test_;
  V2Test v2test_;
  MTest mtest_;
  
  /*! Set up routine.
   */
  virtual void SetUp()
  {
    KernelTestWithType< _ST_ >::SetUp();
    KernelTestWithMap<_N_>::SetUp();
    if (typeImplemented_ && !problemTooSmall_)
    {
      v1test_.SetUp();
      v1test_.replaceMap(map_);
      v2test_.SetUp();
      v2test_.replaceMap(map_);
      // these probably have a different communicator, but this shouldn't matter here...
      mtest_.SetUp();

      V1_    = v1test_.vec1_;
      V1_vp_ = v1test_.vec1_vp_;
      ldaV1_ = v1test_.lda_;

      V2_    = v2test_.vec1_;
      V2_vp_ = v2test_.vec1_vp_;
      ldaV2_ = v2test_.lda_;

      M1_    = mtest_.mat1_;
      M1_vp_ = mtest_.mat1_vp_;
      ldaM1_ = mtest_.m_lda_;

      M2_    = mtest_.mat2_;
      M2_vp_ = mtest_.mat2_vp_;
      ldaM2_ = mtest_.m_lda_;
    }
    else
    {
      V1_ = NULL;
      V1_vp_ = NULL;
      V2_ = NULL;
      V2_vp_ = NULL;
      M1_ = NULL;
      M1_vp_ = NULL;
      M2_ = NULL;
      M2_vp_ = NULL;
    }
    stride_=1;
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      mtest_.TearDown();
      v2test_.TearDown();
      v1test_.TearDown();
    }
    KernelTestWithMap<_N_>::TearDown();
    KernelTestWithType<_ST_>::TearDown();
  }

};

  // check ones(n,m)'*ones(n,k)=n*ones(m,k), and columns with 1, 2, 3...
  TEST_F(CLASSNAME, mvecT_times_mvec) 
  {
    if (typeImplemented_ && !problemTooSmall_)
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
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatEqual(M1_,(ST)nglob_));
      ASSERT_EQ(0,iflag_);

      // fill rows with 1,2,3,4, ...
      for (int ii=0; ii< nloc_; ii++)
      {
        for (int j=0; j<k_; j++)
          V2_vp_[VIDX(ii,j,ldaV2_)] = -st::one()*(_MT_)(j+1);
        for (int i=0; i<m_; i++)
          V1_vp_[VIDX(ii,i,ldaV1_)] = st::one()*(_MT_)((i+1)*1.0l/nglob_);
      }
      SUBR(mvec_to_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_to_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check result
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < k_; j++)
        {
          _MT_ val = -mt::one()*(i+1)*(j+1);
          ASSERT_NEAR(val, st::real(M1_vp_[MIDX(i,j,ldaM1_)]), _N_*20*mt::eps());
          ASSERT_NEAR(mt::zero(), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]), _N_*20*mt::eps());
        }
      }

      // fill columns with row/i and row/(i+1)
      // exploits sum_(k=1)^n 1/(k*(k+1)) = 1 - 1/(n+1)
      gidx_t ilower;     
      phist_map_get_ilower(map_,&ilower,&iflag_);
      for (int ii=0; ii< nloc_; ii++)
      {
        for (int j=0; j<k_; j++)
          V2_vp_[VIDX(ii,j,ldaV2_)] = (_ST_) -st::one()*(_MT_)((j+1)*1.0l/(ilower+ii+1));
        for (int i=0; i<m_; i++)
          V1_vp_[VIDX(ii,i,ldaV1_)] = (_ST_) st::one()*(_MT_)((i+1)*1.0l/(ilower+ii+2));
      }
      SUBR(mvec_to_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_to_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check result
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < k_; j++)
        {
          _MT_ val = -mt::one()*(i+1)*(j+1)*(1.0l - 1.0l/(_N_+1));
          ASSERT_NEAR(val, st::real(M1_vp_[MIDX(i,j,ldaM1_)]), _N_*20*mt::eps());
          ASSERT_NEAR(mt::zero(), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]), _N_*20*mt::eps());
        }
      }


    }
  }

#if ( _M_ == _K_ )
  // check ones(n,m)'*ones(n,m)=n*ones(m,m), and columns with 1, 2, 3...
  TEST_F(CLASSNAME, mvecT_times_mvec_self) 
  {
    if( typeImplemented_ && !problemTooSmall_ && m_ == k_ )
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V1_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      V1Test::PrintVector(*cout,"ones",V1_vp_,nloc_,ldaV1_,stride_,mpi_comm_);
      MTest::PrintSdMat(*cout,"ones'*ones",M1_vp_,ldaM1_,stride_,mpi_comm_);
#endif
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatEqual(M1_,(ST)nglob_));
      ASSERT_EQ(0,iflag_);

      // fill rows with 1,2,3,4, ...
      for (int ii=0; ii< nloc_; ii++)
      {
        for (int j=0; j<k_; j++)
          V1_vp_[VIDX(ii,j,ldaV1_)] = -st::one()*(_MT_)(j+1);
      }
      SUBR(mvec_to_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V1_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check result
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < m_; j++)
        {
          _MT_ val = mt::one()*(i+1)*(j+1)*_N_;
          ASSERT_NEAR(val, st::real(M1_vp_[MIDX(i,j,ldaM1_)]), _N_*10*mt::eps()*_N_);
          ASSERT_NEAR(mt::zero(), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]), _N_*10*mt::eps()*_N_);
        }
      }
    }
  }

  // check V'V gives same as W=V, V'W with random V
  TEST_F(CLASSNAME, mvecT_times_mvec_self_rand)
  {
    if( typeImplemented_ && !problemTooSmall_ && m_ == k_ )
    {
      // fill V with random numbers and set W=V
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),V1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V1_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      MTest::sdMat_parallel_check(M2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatsEqual(M1_,M2_));
      ASSERT_EQ(0,iflag_);

    }
  }
#endif

  // check ones(n,m)'*ones(n,k)=n*ones(m,k)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, mvecT_times_mvec_prec) 
#else
  TEST_F(CLASSNAME, DISABLED_mvecT_times_mvec_prec) 
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(V2_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

//SUBR(sdMat_print)(M1_, &iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatEqual(M1_,(ST)nglob_));
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // fill rows with 1,2,3,4, ...
      for (int ii=0; ii< nloc_; ii++)
      {
        for (int j=0; j<k_; j++)
          V2_vp_[VIDX(ii,j,ldaV2_)] = -st::one()*(_MT_)(j+1);
        for (int i=0; i<m_; i++)
          V1_vp_[VIDX(ii,i,ldaV1_)] = st::one()*(_MT_)((i+1)*1./nglob_);
      }
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check result
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < k_; j++)
        {
          _MT_ val = -mt::one()*(i+1)*(j+1);
          ASSERT_REAL_EQ(val, st::real(M1_vp_[MIDX(i,j,ldaM1_)]));
          ASSERT_REAL_EQ(mt::zero(), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]));
        }
      }

      // fill columns with row/i and row/(i+1)
      // exploits sum_(k=1)^n 1/(k*(k+1)) = 1 - 1/(n+1)
      gidx_t ilower;     
      phist_map_get_ilower(map_,&ilower,&iflag_);
      for (int ii=0; ii< nloc_; ii++)
      {
        for (int j=0; j<k_; j++)
          V2_vp_[VIDX(ii,j,ldaV2_)] = (_ST_) -st::one()*(_MT_)((j+1)*1.0l/(ilower+ii+1));
        for (int i=0; i<m_; i++)
          V1_vp_[VIDX(ii,i,ldaV1_)] = (_ST_) st::one()*(_MT_)((i+1)*1.0l/(ilower+ii+2));
      }
      SUBR(mvec_to_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_to_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check result
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < k_; j++)
        {
          _MT_ val = -mt::one()*(i+1)*(j+1)*(1.0l - 1.0l/(_N_+1));
          ASSERT_NEAR(val, st::real(M1_vp_[MIDX(i,j,ldaM1_)]), 100*mt::eps());
          ASSERT_NEAR(mt::zero(), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]), 100*mt::eps());
        }
      }


    }
  }

#if ( _M_ == _K_ )
  // check ones(n,m)'*ones(n,k)=n*ones(m,k)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, mvecT_times_mvec_prec_self) 
#else
  TEST_F(CLASSNAME, DISABLED_mvecT_times_mvec_prec_self) 
#endif
  {
    if( typeImplemented_ && !problemTooSmall_ && m_ == k_ )
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V1_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

//SUBR(sdMat_print)(M1_, &iflag_);
      ASSERT_REAL_EQ(mt::one(),SdMatEqual(M1_,(ST)nglob_));
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // fill rows with 1,2,3,4, ...
      for (int ii=0; ii< nloc_; ii++)
      {
        for (int j=0; j<k_; j++)
          V1_vp_[VIDX(ii,j,ldaV1_)] = -st::one()*(_MT_)(j+1);
      }
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V1_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check result
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < m_; i++)
      {
        for(int j = 0; j < k_; j++)
        {
          _MT_ val = mt::one()*(i+1)*(j+1)*_N_;
          ASSERT_REAL_EQ(val, st::real(M1_vp_[MIDX(i,j,ldaM1_)]));
          ASSERT_REAL_EQ(mt::zero(), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]));
        }
      }
    }
  }
#endif


  // check ones(n,m)'*ones(n,m)=n*ones(m,m)
  TEST_F(CLASSNAME, mvecT_times_mvec_with_manual_comparison) 
  {
    if (typeImplemented_ && !problemTooSmall_)
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

      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
#ifdef IS_COMPLEX
#define _ST_PREC_ std::complex<long double>
#define CONJ_PREC(x) std::conj(x)
#else
#define _ST_PREC_ long double
#define CONJ_PREC(x) x
#endif
      for (int i=0; i<m_; i++)
      {
        for (int j=0; j<k_; j++)
        {
// melven: Kahan summation here will most probably not work, see remarks below
// * compiler optimization!
// * one can also improve the multiplication
// * high precision examples are in src/kernels/builtin/\*prec\*
// -> long double for reference data should be the better way (hoping/checking that long double != double!)
          _ST_PREC_ dot_ij = (_ST_PREC_)0;
          long double dotAbs_ij = (long double)0;
          for (int ii=0; ii< nloc_; ii++)
          {
            _ST_PREC_ tmp = CONJ_PREC(V1_vp_[VIDX(ii,i,ldaV1_)])*V2_vp_[VIDX(ii,j,ldaV2_)];
            dot_ij += tmp;
            dotAbs_ij += std::abs(tmp);
          }
#ifdef PHIST_HAVE_MPI
#ifdef IS_COMPLEX
          MPI_Allreduce(MPI_IN_PLACE,&dot_ij,2,MPI_LONG_DOUBLE,MPI_SUM,mpi_comm_);
#else
          MPI_Allreduce(MPI_IN_PLACE,&dot_ij,1,MPI_LONG_DOUBLE,MPI_SUM,mpi_comm_);
#endif
          MPI_Allreduce(MPI_IN_PLACE,&dotAbs_ij,1,MPI_LONG_DOUBLE,MPI_SUM,mpi_comm_);
#endif
          _MT_ cond = dotAbs_ij / std::abs(dot_ij);
          EXPECT_NEAR(st::real((_ST_)dot_ij), st::real(M1_vp_[MIDX(i,j,ldaM1_)]), dotAbs_ij*cond*1000*mt::eps());
          EXPECT_NEAR(st::imag((_ST_)dot_ij), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]), dotAbs_ij*cond*1000*mt::eps());
        }
      }
    }
#undef _ST_PREC_
#undef CONJ_PREC
  }


  // check ones(n,m)'*ones(n,m)=n*ones(m,m)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, mvecT_times_mvec_with_manual_comparison_prec_hard)
#else
  TEST_F(CLASSNAME, DISABLED_mvecT_times_mvec_with_manual_comparison_prec_hard)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with random numbers
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_from_device)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
#ifdef IS_COMPLEX
#define _ST_PREC_ std::complex<long double>
#define CONJ_PREC(x) std::conj(x)
#else
#define _ST_PREC_ long double
#define CONJ_PREC(x) x
#endif
      for (int i=0; i<m_; i++)
      {
        for (int j=0; j<k_; j++)
        {
// melven: Kahan summation here will most probably not work, see remarks below
// * compiler optimization!
// * one can also improve the multiplication
// * high precision examples are in src/kernels/builtin/*prec*
// -> long double for reference data should be the better way (hoping/checking that long double != double!)
          _ST_PREC_ dot_ij = (_ST_PREC_)0;
          long double dotAbs_ij = (long double)0;
          for (int ii=0; ii< nloc_; ii++)
          {
            _ST_PREC_ tmp = CONJ_PREC(V1_vp_[VIDX(ii,i,ldaV1_)])*V2_vp_[VIDX(ii,j,ldaV2_)];
            dot_ij += tmp;
            dotAbs_ij += std::abs(tmp);
          }
#ifdef PHIST_HAVE_MPI
#ifdef IS_COMPLEX
          MPI_Allreduce(MPI_IN_PLACE,&dot_ij,2,MPI_LONG_DOUBLE,MPI_SUM,mpi_comm_);
#else
          MPI_Allreduce(MPI_IN_PLACE,&dot_ij,1,MPI_LONG_DOUBLE,MPI_SUM,mpi_comm_);
#endif
          MPI_Allreduce(MPI_IN_PLACE,&dotAbs_ij,1,MPI_LONG_DOUBLE,MPI_SUM,mpi_comm_);
#endif
          _MT_ cond = dotAbs_ij / std::abs(dot_ij);
          PHIST_SOUT(PHIST_INFO, "error: %e (cond. number: %e, eps: %e)\n", st::abs((_ST_)dot_ij-M1_vp_[MIDX(i,j,ldaM1_)]),cond,mt::eps());
          EXPECT_NEAR(st::real((_ST_)dot_ij), st::real(M1_vp_[MIDX(i,j,ldaM1_)]), 20*cond*mt::eps());
          EXPECT_NEAR(st::imag((_ST_)dot_ij), st::imag(M1_vp_[MIDX(i,j,ldaM1_)]), 20*cond*mt::eps());
        }
      }
    }
#undef _ST_PREC_
#undef CONJ_PREC
  }

  // check ones(n,m)*ones(m,k)=m*ones(n,k),
  // and ones(n,m)*ones(m,k)-m*ones(n,k)=0
  TEST_F(CLASSNAME, mvec_times_sdMat)
  {
    if (typeImplemented_ && !problemTooSmall_)
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
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // with adding to factor*(input vector)
      ST alpha=(ST)2.0+(ST)2.0*st::cmplx_I();
      SUBR(mvec_put_value)(V2_,(ST)m_*alpha,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,-st::one()/alpha,V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V2_,st::zero()));

      // with adding factor*V*C to input vector
      SUBR(mvec_put_value)(V2_,(MT)m_*alpha,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(-alpha,V1_,M1_,st::one(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V2_,st::zero()));
    }
  }

  // check ones(n,m)*ones(m,k)=m*ones(n,k),
  // and ones(n,m)*ones(m,k)-m*ones(n,k)=0
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, mvec_times_sdMat_prec)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_prec)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(V2_,(ST)42.0*st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
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
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // with adding to factor*(input vector)
      ST alpha=(ST)2+(ST)2*st::cmplx_I();
      SUBR(mvec_put_value)(V2_,(ST)m_*alpha,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,-st::one()/alpha,V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V2_,st::zero()));

      // with adding factor*V*C to input vector
      SUBR(mvec_put_value)(V2_,(MT)m_*alpha,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(-alpha,V1_,M1_,st::one(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V2_,st::zero()));
    }
  }

#if (_M_>=_K_)
  // check ones(n,m)*ones(m,m)=m*ones(n,m)
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place)
  {
    if (typeImplemented_ && !problemTooSmall_)
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
      TYPE(mvec_ptr) V_cols=NULL;
      SUBR(mvec_view_block)(V1_,&V_cols,0,k_-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(V_cols,(ST)m_));
   
      if (k_<m_)
      {
        SUBR(mvec_view_block)(V1_,&V_cols,k_,m_-1,&iflag_);
      }
      SUBR(mvec_delete)(V_cols,&iflag_);
      ASSERT_EQ(0,iflag_);

      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check ones(n,m)*ones(m,m)=m*ones(n,m)
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place_with_random_data)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with ones
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(M1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
      MTest::sdMat_parallel_check(M1_,&iflag_);
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
  
      TYPE(mvec_ptr) V_cols=NULL;
      SUBR(mvec_view_block)(V1_,&V_cols,0,k_-1,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MvecsEqual(V_cols,V2_),sqrt(mt::eps()));
   
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_delete)(V_cols,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check ones(n,m)*ones(m,m)=m*ones(n,m)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place_prec)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_in_place_prec)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_ && m_ == k_ )
    {
      // fill V and W with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_put_value)(M1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(V1_,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) V_cols=NULL;
      SUBR(mvec_view_block)(V1_,&V_cols,0,k_-1,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),MvecEqual(V_cols,(ST)m_));
   
      SUBR(mvec_delete)(V_cols,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // check ones(n,m)*ones(m,m)=m*ones(n,m)
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place_with_random_data_prec)
#else
  TEST_F(CLASSNAME, DISABLED_mvec_times_sdMat_in_place_with_random_data_prec)
#endif
  {
    if (typeImplemented_ && !problemTooSmall_ && m_ == k_ )
    {
      // fill V and W with ones
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_sync_values)(M1_, comm_, &iflag_);
      ASSERT_EQ(0,iflag_);
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat)(st::one(),V1_,M1_,st::zero(),V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_times_sdMat_inplace)(V1_,M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
  
      TYPE(mvec_ptr) V_cols=NULL;
      SUBR(mvec_view_block)(V1_,&V_cols,0,k_-1,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MvecsEqual(V_cols,V2_),sqrt(mt::eps()));
   
      SUBR(mvec_delete)(V_cols,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


#if (_M_==_K_)
  // check that we can zero out some columns of V by multiplying
  // a view of them with a zero sdMat.
  TEST_F(CLASSNAME, mvec_times_sdMat_in_place_with_views)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and M with ones
      SUBR(mvec_put_value)(V1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(V1_,&iflag_);
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

      SUBR(sdMat_delete)(M,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(V,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

/* _M_==_K_ */
#endif

/* _M_>=_K_ */
#endif

  // random check
  TEST_F(CLASSNAME, random_mvecT_times_mvec) 
  {
    if (typeImplemented_ && !problemTooSmall_)
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
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
      // to check the result, scale V1 by -1, add V1'V2 again and compare with 0
      SUBR(mvec_scale)(V1_,-st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::one(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
  
      // but for large cases we need to respect the condition number as it won't be zero at all
      _MT_ normV1[_M_];
      SUBR(mvec_norm2)(V1_, normV1, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ normV2[_K_];
      SUBR(mvec_norm2)(V2_, normV2, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ maxNorm12 = mt::zero();
      for(int i = 0; i < _M_; i++)
        for(int j = 0; j < _K_; j++)
          maxNorm12 = std::max(maxNorm12, normV1[i]*normV2[j]);

      ASSERT_NEAR(mt::one(),SdMatEqual(M1_,st::zero()),maxNorm12*100*mt::eps());
    
    }
  }

  // random check
#if defined(PHIST_HIGH_PRECISION_KERNELS) && !USE_VIEWS
  TEST_F(CLASSNAME, random_mvecT_times_mvec_prec) 
#else
  TEST_F(CLASSNAME, DISABLED_random_mvecT_times_mvec_prec) 
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      // fill V and W with ones
      SUBR(mvec_random)(V1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(V2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::zero(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    
      // to check the result, scale V1 by -1, add V1'V2 again and compare with 0
      SUBR(mvec_scale)(V1_,-st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);

      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvecT_times_mvec)(st::one(),V1_,V2_,st::one(),M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
  
      // but for large cases we need to respect the condition number as it won't be zero at all
      _MT_ normV1[_M_];
      SUBR(mvec_norm2)(V1_, normV1, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ normV2[_K_];
      SUBR(mvec_norm2)(V2_, normV2, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ maxNorm12 = mt::zero();
      for(int i = 0; i < _M_; i++)
        for(int j = 0; j < _K_; j++)
          maxNorm12 = std::max(maxNorm12, normV1[i]*normV2[j]);

      // TODO: what is the correctly required precision here? 
      ASSERT_NEAR(mt::one(),SdMatEqual(M1_,st::zero()),mt::sqrt(maxNorm12)*100*mt::eps());
    }
  }

  // random check with partial views of partial mvecs and sdMats
  TEST_F(CLASSNAME, random_mvecT_times_mvec_with_inside_views)
  {
    if (typeImplemented_ && !problemTooSmall_ && m_ > 4)
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

      // for large cases we need to respect the condition number (which requires sum |v1_i*v2_i| bounded by the ||v_1|| * ||v_2||)
      _MT_ normV1[_M_];
      SUBR(mvec_norm2)(V1_, normV1, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ normV2[_K_];
      SUBR(mvec_norm2)(V2_, normV2, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ maxNorm12 = mt::one();
      for(int i = 0; i < _M_; i++)
        for(int j = 0; j < _K_; j++)
          maxNorm12 = std::max(maxNorm12, normV1[i]*normV2[j]);

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
        MTest::sdMat_parallel_check(M1_,&iflag_);
        ASSERT_EQ(0,iflag_);

        MTest::sdMat_parallel_check(M2_,&iflag_);
        ASSERT_EQ(0,iflag_);

        MTest::sdMat_parallel_check(M3,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        // subtract the result without views from the one with non-viewed target sdMat
        SUBR(sdMat_add_sdMat)(-st::one(),M1, st::one(),M3,&iflag_);
        ASSERT_EQ(0,iflag_);

        // the result should be zero!
        SUBR(sdMat_from_device)(M3,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),ArrayEqual(M3_vp,m1[i],m2[i],lda_M3,stride_,st::zero(),mflag_),_N_*maxNorm12*10*mt::eps());

        // subtract the two viewed blocks in the result sdMats
        SUBR(sdMat_add_sdMat)(-st::one(),M1, st::one(),M2,&iflag_);
        ASSERT_EQ(0,iflag_);

        // the result should be zero!
        SUBR(sdMat_from_device)(M2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),ArrayEqual(M2_vp_,m_,k_,ldaM2_,stride_,st::zero(),mflag_),_N_*maxNorm12*10*mt::eps());

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

#ifdef PHIST_KERNEL_LIB_BUILTIN
  // the result of sdMat_random must be equal on all processes even after calling mvec_random
  TEST_F(CLASSNAME, parallel_random)
  {
    if (typeImplemented_ && !problemTooSmall_ )
    {
      SUBR(mvec_random)(V1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // and do it again
      SUBR(mvec_random)(V1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_random)(M1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      MTest::sdMat_parallel_check(M1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
#endif

const int CLASSNAME::k_;
const int CLASSNAME::m_;
const int CLASSNAME::n_;
