#include "../tools/TestHelpers.h"
#include <algorithm>

#ifndef CLASSNAME
#error "file not included correctly"
#endif

#ifdef F_INIT
#undef F_INIT
#endif
//#define F_INIT(a,b)  (ONE/((_ST_)(a)+(_ST_)(b)+ONE) + _CMPLX_I_*((_ST_)(a)-(_ST_)(b)))
#define F_INIT(a,b)  (st::one()/((ST)(a)+(ST)(b)+st::one()) +st::cmplx_I()*((ST)(a)-(ST)(b)))


// some function to test mvec_put_func
#ifdef FIRST_INSTANCE
int PREFIX(mvecInitializer)(ghost_gidx_t i, ghost_lidx_t j, void* vval)
{
#include "phist_std_typedefs.hpp"
  _ST_* val = (_ST_*)vval;
  *val = F_INIT(i,j);
  return 0;
}
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

      // we need to somehow initialize the vectors. Here we choose
      // to do it manually (and not rely on the functions mvec_random
      // and mvec_put_value provided by the kernel lib).
      // In order for the tests to work we must get
      // the data to the GPU if necessary, even if this
      // functionality has not been tested (strictly speaking).
      // In SdMat_test we use the provided functions instead, hopefully
      // the diversity here will lead to a well-covered code in total.
            
      for (int j=0;j<nvec_;j++)
        for (int i=0;i<nloc_*stride_;i+=stride_)
          {
          vec1_vp_[VIDX(i,j,lda_)]=random_number();
          vec2_vp_[VIDX(i,j,lda_)]=st::one();
          }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      SUBR(mvec_to_device)(vec2_,&iflag_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    }

};

  TEST_F(CLASSNAME, my_length)
    {
    if (typeImplemented_)
      {
      lidx_t nloc;
      SUBR(mvec_my_length)(vec1_,&nloc,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nloc_, nloc); 
      }
    }

  TEST_F(CLASSNAME, num_vectors) 
    {
    if (typeImplemented_)
      {
      int nvec;
      SUBR(mvec_num_vectors)(vec1_,&nvec,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(nvec_, nvec);
      }
    }


  TEST_F(CLASSNAME, put_value) 
    {
    if (typeImplemented_)
      {
      ST val = (_ST_)42.0 + (ST)3.0*st::cmplx_I();
      SUBR(mvec_put_value)(vec1_,val,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec1_,val));

      // check that the put_value function does not change the pointer
      ST* ptr;
      lidx_t lda_new;
      SUBR(mvec_extract_view)(vec1_,&ptr,&lda_new,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda_,lda_new);
      ASSERT_EQ(vec1_vp_,ptr);
      }
    }

  TEST_F(CLASSNAME, dot_mvec)
  {
    if (typeImplemented_)
    {
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
        for (int i=0;i<nloc_*stride_;i+=stride_)
        {
          vec2_vp_[VIDX(i,j,lda_)]=mt::one()/st::conj(vec1_vp_[VIDX(i,j,lda_)]);
        }
      SUBR(mvec_to_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ dots_ref[_NV_];
      _ST_ dots[_NV_];
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots_ref,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ val = st::one() * (ST)nglob_;
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(dots_ref,nvec_,1,nvec_,1,val));

      // test two random vectors
      SUBR(mvec_random)(vec2_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots_ref,&iflag_);
      ASSERT_EQ(0,iflag_);

      _MT_ dotsAbs[_NV_];
      for (int j=0;j<nvec_;j++)
      {
        dots[j] = st::zero();
        dotsAbs[j] = mt::zero();
        for (int i=0;i<nloc_*stride_;i+=stride_)
        {
          _ST_ tmp = st::conj(vec1_vp_[VIDX(i,j,lda_)])*vec2_vp_[VIDX(i,j,lda_)];
          dots[j] += tmp;
          dotsAbs[j] += st::abs(tmp);
        }
#ifdef PHIST_HAVE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &dots[j], 1, st::mpi_type(), MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &dotsAbs[j], 1, mt::mpi_type(), MPI_SUM, MPI_COMM_WORLD);
#endif
        _MT_ cond = dotsAbs[j] / st::abs(dots[j]);
        EXPECT_NEAR(mt::zero(), st::real(dots[j]-dots_ref[j]), dotsAbs[j]*10*mt::eps());
        EXPECT_NEAR(mt::zero(), st::imag(dots[j]-dots_ref[j]), dotsAbs[j]*10*mt::eps());
      }
    }

  }

#if (_N_ % 4 == 0)
  // mvec_dot_mvec tests, more precise version
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, dot_mvec_prec)
#else
  TEST_F(CLASSNAME, DISABLED_dot_mvec_prec)
#endif
  {
    if (typeImplemented_)
    {
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
        for (int i=0;i<nloc_*stride_;i+=stride_)
        {
          vec2_vp_[VIDX(i,j,lda_)]=mt::one()/st::conj(vec1_vp_[VIDX(i,j,lda_)]);
        }
      SUBR(mvec_to_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ dots_ref[_NV_];
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots_ref,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ val = st::one() * (ST)nglob_;
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(dots_ref,nvec_,1,nvec_,1,val));

      // test two random vectors
      SUBR(mvec_random)(vec2_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots_ref,&iflag_);
      ASSERT_EQ(0,iflag_);

#ifdef IS_COMPLEX
#define _ST_PREC_ std::complex<long double>
#define CONJ_PREC(x) std::conj(x)
#else
#define _ST_PREC_ long double
#define CONJ_PREC(x) x
#endif
      long double dotsAbs[_NV_];
      _ST_PREC_ dots[_NV_];
      for (int j=0;j<nvec_;j++)
      {
        dots[j] = (_ST_PREC_)0;
        dotsAbs[j] = (long double)0;
        for (int i=0;i<nloc_*stride_;i+=stride_)
        {
          _ST_PREC_ tmp = CONJ_PREC(vec1_vp_[VIDX(i,j,lda_)])*vec2_vp_[VIDX(i,j,lda_)];
          dots[j] += tmp;
          dotsAbs[j] += std::abs(tmp);
        }
#ifdef PHIST_HAVE_MPI
#ifdef IS_COMPLEX
        MPI_Allreduce(MPI_IN_PLACE, &dots[j], 2, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        MPI_Allreduce(MPI_IN_PLACE, &dots[j], 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        MPI_Allreduce(MPI_IN_PLACE, &dotsAbs[j], 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        // cond. number of the summation
        _MT_ cond = (_MT_)(dotsAbs[j] / std::abs(dots[j]));
        PHIST_SOUT(PHIST_INFO, "error: %e (cond. number: %e, eps: %e)\n", st::abs((_ST_)dots[j]-dots_ref[j]),cond,mt::eps());
        EXPECT_NEAR(mt::zero(), st::real((_ST_)dots[j]-dots_ref[j]), 2*cond*mt::eps());
        EXPECT_NEAR(mt::zero(), st::imag((_ST_)dots[j]-dots_ref[j]), 2*cond*mt::eps());
      }
    }
#undef _ST_PREC_
#undef CONJ_PREC

  }

  // mvec_dot_mvec tests, more precise version
  // exploits sum_(k=1)^n 1/(k*(k+1)) = 1 - 1/(n+1)
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, dot_mvec_prec_hard)
#else
  TEST_F(CLASSNAME, DISABLED_dot_mvec_prec_hard)
#endif
  {
    if (typeImplemented_)
    {
      gidx_t ilower;     
      phist_map_get_ilower(map_,&ilower,&iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
      {
        for(int i = 0; i < nloc_; i++)
        {
          vec2_vp_[VIDX(i,j,lda_)]=st::one()*((MT)1./(ilower+i+1));
          vec1_vp_[VIDX(i,j,lda_)]=st::one()*((MT)1./(ilower+i+2));
        }
      }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_to_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ dots[_NV_];
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ val = (_ST_)(st::one()*(_MT_)(1.0l-1.0l/(_N_+1)));
      PHIST_SOUT(PHIST_INFO, "error: %e\n", st::abs(dots[0]-val));
      // melven: this is a *hard* precision test
      // DO NOT CHANGE THIS TO ASSERT_NEAR!
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(dots,nvec_,1,nvec_,1,val));
    }
  }

  // mvec_dot_mvec tests, more precise version
  // half of the numbers a 1, then one 1.e-18, then the other half is -1
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, dot_mvec_prec_veryhard)
#else
  TEST_F(CLASSNAME, DISABLED_dot_mvec_prec_veryhard)
#endif
  {
    if (typeImplemented_)
    {
      gidx_t ilower;     
      phist_map_get_ilower(map_,&ilower,&iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
      {
        for(int i = 0; i < nloc_; i++)
        {
          if( ilower+i+1 < _N_/2 )
            vec2_vp_[VIDX(i,j,lda_)]=st::one();
          else if( ilower+i+1 == _N_/2 || ilower+i == _N_/2)
            vec2_vp_[VIDX(i,j,lda_)]=st::eps();
          else
            vec2_vp_[VIDX(i,j,lda_)]=-st::one();
          vec1_vp_[VIDX(i,j,lda_)]=st::one();
        }
      }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_to_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ dots[_NV_];
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ val = 2*st::eps();
      PHIST_SOUT(PHIST_INFO, "error: %e\n", st::abs(dots[0]-val));
      // melven: this is a *hard* precision test
      // DO NOT CHANGE THIS TO ASSERT_NEAR!
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(dots,nvec_,1,nvec_,1,val));
    }
  }
#endif


#if _N_ < 100
  TEST_F(CLASSNAME, random)
    {
    if (typeImplemented_)
      {
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check that the random function does not change the pointer
      ST* ptr;
      lidx_t lda_new;
      SUBR(mvec_extract_view)(vec1_,&ptr,&lda_new,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(lda_,lda_new);
      ASSERT_EQ(vec1_vp_,ptr);
            
      MT absval[nloc_*nvec_];
      int k=0;
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<nloc_*stride_;i+=stride_)
          {
          absval[k++]=st::abs(vec1_vp_[VIDX(i,j,lda_)]);
          }
        }
      MT minval=1.0;
      std::sort(absval,absval+k);
      for (int j=1;j<k;j++)
        {
        minval=std::min(minval,mt::abs(absval[j]-absval[j-1]));
        }
      // force assertion failure if two 'random' numbers are the same
      // TODO - this test sometimes fails with OMP_NUM_THREADS>1 and ghost
      //         in single precision. Probably the chance to get 'identical'
      //         numbers in SP is higher, but why does it depend on #threads?
      if (minval<=mt::eps())
      {
        SUBR(mvec_print)(vec1_,&iflag_);
      }
      ASSERT_EQ(true,minval>mt::eps()); 
      }
    }
#endif

  TEST_F(CLASSNAME, upload_download)
  {
    // just tests that the upload and from_device functions return 0
    // and do not crash.
    if (typeImplemented_)
    {
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  // 2-norm, nrm2=sqrt(v'v)
  TEST_F(CLASSNAME, norm2)
  {
    if (typeImplemented_)
    {
      gidx_t ilower;     
      phist_map_get_ilower(map_,&ilower,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
      {
        for (int i=0;i<nloc_*stride_;i+=stride_)
        {
          vec1_vp_[VIDX(i,j,lda_)]=ilower+i;
        }
      }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // sum_(k=0)^(n-1) k^2 = 1/6*(n-1)*n*(2*n-1)
      MT expect = (MT)std::sqrt(1.0l/6.0l*(_N_-1.0l)*_N_*(2.0l*_N_-1.0l));

      MT nrm2[nvec_];
      SUBR(mvec_norm2)(vec1_,nrm2,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int i=0;i<nvec_;i++)
      {
        ASSERT_NEAR(expect,nrm2[i],expect*mt::sqrt(mt::eps()));
      }
    }
  }

#if (_N_ % 4 == 0 )
  // 2-norm, nrm2=sqrt(v'v), more precise version
#ifdef PHIST_HIGH_PRECISION_KERNELS
  TEST_F(CLASSNAME, norm2_precise)
#else
  TEST_F(CLASSNAME, DISABLED_norm2_precise)
#endif
  {
    if (typeImplemented_)
    {
      gidx_t ilower;     
      phist_map_get_ilower(map_,&ilower,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
      {
        for (int i=0;i<nloc_*stride_;i+=stride_)
        {
          vec1_vp_[VIDX(i,j,lda_)]=ilower+i;
        }
      }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // sum_(k=0)^(n-1) k^2 = 1/6*(n-1)*n*(2*n-1)
      MT expect = (MT)std::sqrt(1.0l/6.0l*(_N_-1.0l)*_N_*(2.0l*_N_-1.0l));

      MT nrm2[nvec_];
      iflag_ = PHIST_ROBUST_REDUCTIONS;
      SUBR(mvec_norm2)(vec1_,nrm2,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int i=0;i<nvec_;i++)
      {
        ASSERT_REAL_EQ(expect,nrm2[i]);
      }
    }
  }
#endif

  // 2-norm, nrm2=sqrt(v'v)
  TEST_F(CLASSNAME, norm2_of_viewed_cols)
  {
    if (typeImplemented_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      MT nrm2[nvec_];
      SUBR(mvec_norm2)(vec1_,nrm2,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int jmin=std::min(1,nvec_-1);
      int jmax=std::min(3,nvec_-1);
      TYPE(mvec_ptr) view=NULL;
      SUBR(mvec_view_block)(vec1_,&view,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);
      MT nrm2_view[jmax-jmin+1];

      SUBR(mvec_norm2)(view,nrm2_view,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      for (int j=jmin;j<=jmax;j++)
      {
        ASSERT_NEAR(nrm2[j],nrm2_view[j-jmin],_N_*10*mt::eps());
      }
      
      SUBR(mvec_delete)(view,&iflag_);
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
      SUBR(mvec_add_mvec)(alpha,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec1_,vec2_));
      }
    
    }

  // X = 0*Y + a*X = a*X
  TEST_F(CLASSNAME, scale_by_axpy)
    {
    if (typeImplemented_)
      {
      ST alpha = st::zero();
      ST beta  = st::prand();
      PHIST_OUT(9,"axpy, alpha=%f+%f i, beta=%f+%f i",st::real(alpha),
        st::imag(alpha),st::real(beta),st::imag(beta));
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(vec2_,&iflag_);
      PrintVector(std::cerr,"before scale",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
#endif
      SUBR(mvec_add_mvec)(alpha,vec1_,beta,vec2_,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_from_device)(vec2_,&iflag_);
      PrintVector(std::cerr,"after scale",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
#endif
      ASSERT_EQ(0,iflag_);
            
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec2_,beta));
      }
    
    }

  // X = Y*diag(a_1,...,a_nvec) + a*X
  TEST_F(CLASSNAME, random_add)
  {
    if( typeImplemented_ )
    {
      ST beta = st::prand();
      ST alpha= st::prand();

      SUBR(mvec_add_mvec)(alpha,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // calculate solution by hand
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          vec1_vp_[VIDX(i,j,lda_)] = alpha*vec1_vp_[VIDX(i,j,lda_)]+beta;
        }
      }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),MvecsEqual(vec1_,vec2_), 1000*mt::eps());
    }
  }

  // X = Y*diag(a_1,...,a_nvec) + a*X
  TEST_F(CLASSNAME, random_vadd)
  {
    if( typeImplemented_ )
    {
      ST beta = st::prand();
      ST alpha[_NV_];
      for(int i = 0; i < nvec_; i++)
        alpha[i] = st::prand();

      SUBR(mvec_vadd_mvec)(alpha,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // calculate solution by hand
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          vec1_vp_[VIDX(i,j,lda_)] = alpha[j]*vec1_vp_[VIDX(i,j,lda_)]+beta;
        }
      }
      SUBR(mvec_to_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),MvecsEqual(vec1_,vec2_), 1000*mt::eps());
    }
  }

  // view certain columns, manipulate them and check it changes the
  // correct locations in the original one
  TEST_F(CLASSNAME, view_block)
    {
    if (typeImplemented_)
      {
      int jmin=std::min(2,nvec_-1);
      int jmax=std::min(5,nvec_-1);
      TYPE(mvec_ptr) v1_view=NULL;
      SUBR(mvec_view_block)(vec1_,&v1_view,jmin,jmax,&iflag_);
      // create a view of the view
      TYPE(mvec_ptr) v1_vv=NULL;
      SUBR(mvec_view_block)(v1_view,&v1_vv,0,jmax-jmin,&iflag_);
      ASSERT_EQ(0,iflag_);

      // now this should delete the original view and create a new one,
      // all vectors must remain valid:
      SUBR(mvec_view_block)(vec1_,&v1_view,jmin,jmax,&iflag_);
      
      _MT_ norms_V1[nvec_];
      _MT_ norms_V1view[nvec_];
      _MT_ norms_V1vv[nvec_];
      
      SUBR(mvec_norm2)(vec1_,norms_V1,&iflag_);
      ASSERT_EQ(0,iflag_);      

      SUBR(mvec_norm2)(v1_view,norms_V1view,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_norm2)(v1_vv,norms_V1vv,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // compare elements one-by-one because our ArraysEqual expects ST rather than MT as 
      // type here.
      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_NEAR(norms_V1[j],norms_V1view[j-jmin],_N_*10*mt::eps());
        ASSERT_NEAR(norms_V1[j],norms_V1vv[j-jmin],_N_*10*mt::eps());
        }
      // set all the viewed entries to a certain value and check that the original vector is 
      // changed.
      _ST_ val = random_number();
      SUBR(mvec_put_value)(v1_view,val,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(mt::one(),ArrayEqual(&vec1_vp_[VIDX(0,jmin,lda_)],nloc_,jmax-jmin+1,lda_,stride_,val,vflag_),100*mt::eps());

      // new norms after changing columns
      SUBR(mvec_norm2)(vec1_,norms_V1,&iflag_);
      ASSERT_EQ(0,iflag_);      
      
      // delete the views without harming v1
      SUBR(mvec_delete)(v1_view,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(v1_vv,&iflag_);
      ASSERT_EQ(0,iflag_);
      // check that v1 still works
      SUBR(mvec_norm2)(vec1_,norms_V1vv,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int j=0;j<nvec_;j++)
        {
        ASSERT_NEAR(norms_V1[j],norms_V1vv[j],_N_*10*mt::eps());
        }
      }
    }


  // create a view of a view and check if this behaves as the user would suspect
  // (not knowing wether a mvec is actually a view or not!)
  TEST_F(CLASSNAME, nested_view_block)
  {
    if (typeImplemented_)
    {
      // first set some data of the whole array
      _ST_ outer_val = st::prand();
      SUBR(mvec_put_value)(vec1_,outer_val,&iflag_);
      ASSERT_EQ(0,iflag_);

      // now create a view
      int jmin=std::min(2,nvec_-1);
      int jmax=std::min(5,nvec_-1);
      TYPE(mvec_ptr) view = NULL;
      SUBR(mvec_view_block)(vec1_,&view,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      // set the data in the view to some other value
      _ST_ view_val = st::prand();
      SUBR(mvec_put_value)(view,view_val,&iflag_);
      ASSERT_EQ(0,iflag_);

      // view part of view
      int nvec_view;
      SUBR(mvec_num_vectors)(view, &nvec_view, &iflag_);
      ASSERT_EQ(0,iflag_);
      int jmin2=std::min(1,nvec_view-1);
      int jmax2=std::min(3,nvec_view-1);
      TYPE(mvec_ptr) view2 = NULL;
      SUBR(mvec_view_block)(view, &view2, jmin2, jmax2, &iflag_);
      ASSERT_EQ(0,iflag_);

      // set data in the inner view to yet another value
      _ST_ inner_val = st::prand();
      SUBR(mvec_put_value)(view2, inner_val, &iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // now use the raw data to verify results
      _MT_ outerErr = mt::zero();
      _MT_ viewErr = mt::zero();
      _MT_ innerErr = mt::zero();
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          if( j < jmin || j > jmax )
          {
            outerErr = std::max(outerErr, st::abs(vec1_vp_[VIDX(i,j,lda_)]-outer_val));
          }
          else if( j < jmin+jmin2 || j > jmin+jmax2 )
          {
            viewErr = std::max(viewErr, st::abs(vec1_vp_[VIDX(i,j,lda_)]-view_val));
          }
          else
          {
            innerErr = std::max(innerErr, st::abs(vec1_vp_[VIDX(i,j,lda_)]-inner_val));
          }
        }
      }
      ASSERT_REAL_EQ(mt::zero(), outerErr);
      ASSERT_REAL_EQ(mt::zero(), viewErr);
      ASSERT_REAL_EQ(mt::zero(), innerErr);

      SUBR(mvec_delete)(view2, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(view, &iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  // copy in and out columns
  TEST_F(CLASSNAME, get_set_block)
    {
    if (typeImplemented_)
      {
      int jmin=std::min(2,nvec_-1);
      int jmax=std::min(5,nvec_-1);
      TYPE(mvec_ptr) v1_copy=NULL;
      PHISTTEST_MVEC_CREATE(&v1_copy,map_,jmax-jmin+1,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(mvec_get_block)(vec1_,v1_copy,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      _MT_ norms_V1[nvec_];
      _MT_ norms_V1copy[nvec_];
      SUBR(mvec_norm2)(vec1_,norms_V1,&iflag_);
      ASSERT_EQ(0,iflag_);      

      SUBR(mvec_norm2)(v1_copy,norms_V1copy,&iflag_);
      ASSERT_EQ(0,iflag_);

      // compare elements one-by-one because our ArraysEqual expects ST rather than MT as 
      // type here.
      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_NEAR(norms_V1[j],norms_V1copy[j-jmin],_N_*10*mt::eps());
        }
      // set all the viewed entries to a certain value and check that the original vector is 
      // changed.
      _ST_ val = random_number();
      SUBR(mvec_put_value)(v1_copy,val,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check that the norms of v1 are unchanged
      SUBR(mvec_norm2)(vec1_,norms_V1copy,&iflag_);
      ASSERT_EQ(0,iflag_);

      for (int j=0;j<nvec_;j++)
        {
        ASSERT_NEAR(norms_V1[j],norms_V1copy[j],_N_*10*mt::eps());
        }

      // compute the new norms
      // check that the norms of v1 are unchanged
      SUBR(mvec_norm2)(v1_copy,norms_V1copy,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // now set the block as the corresponding columns of v2
      SUBR(mvec_set_block)(vec2_,v1_copy,jmin,jmax,&iflag_);
      ASSERT_EQ(0,iflag_);

      // compute the new norms
      // check that the norms of v1 are unchanged
      SUBR(mvec_norm2)(vec2_,norms_V1,&iflag_);
      ASSERT_EQ(0,iflag_);

      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_NEAR(norms_V1[j],norms_V1copy[j-jmin], _N_*10*mt::eps());
        }

      SUBR(mvec_delete)(v1_copy,&iflag_);
      ASSERT_EQ(0,iflag_);
      }
    }

  TEST_F(CLASSNAME, scale)
  {
    if( typeImplemented_ )
    {
      _ST_ scale = st::prand();

      SUBR(mvec_scale)(vec2_,scale,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec2_,scale));

      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      scale = st::prand();
      SUBR(mvec_scale)(vec1_,scale,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // apply scale to vec2_ by hand
      for(int i = 0; i < nloc_; i++)
        for(int j = 0; j < nvec_; j++)
        {
          vec2_vp_[VIDX(i,j,lda_)] *= scale;
        }
      SUBR(mvec_to_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec1_,vec2_));
    }
  }

  TEST_F(CLASSNAME, vscale)
  {
    if( typeImplemented_ )
    {
      _ST_ scale[_NV_];
      for(int i = 0; i < _NV_; i++)
        scale[i] = st::prand();

      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_vscale)(vec1_,scale,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_from_device)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // apply scale to vec2_ by hand
      for(int i = 0; i < nloc_; i++)
        for(int j = 0; j < nvec_; j++)
        {
          vec2_vp_[VIDX(i,j,lda_)] *= scale[j];
        }

      ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec1_,vec2_));
    }
  }


TEST_F(CLASSNAME,put_func)
{
  if (!typeImplemented_) return;
  SUBR(mvec_put_func)(vec1_,&PREFIX(mvecInitializer),&iflag_);
  ASSERT_EQ(0,iflag_);
  
  gidx_t ilower, iupper;
  phist_map_get_ilower(map_,&ilower,&iflag_);
  ASSERT_EQ(0,iflag_);
  phist_map_get_iupper(map_,&iupper,&iflag_);
  ASSERT_EQ(0,iflag_);
  
  ASSERT_EQ(nloc_,iupper-ilower+1);
  
  for (lidx_t i=0;i<nloc_; i++)
  {
    gidx_t ii = (gidx_t)i + ilower;
    for (lidx_t j=0; j<nvec_; j++)
    {
      vec2_vp_[VIDX(i,j,lda_)]=F_INIT(ii,j);
    }
  }
  SUBR(mvec_to_device)(vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec1_,vec2_));
}

// only test the Belos interface for ghost, we didn't write
// the interfaces for Epetra or Tpetra so it is not our problem.
#ifdef DO_BELOS_TESTS
  // runs all tests from the Belos MvTraits tester
  TEST_F(CLASSNAME, belos_iface)
  {
    if (typeImplemented_)
    {
      Teuchos::RCP<Belos::OutputManager<ST> > MyOM
        = Teuchos::rcp( new Belos::OutputManager<ST>() );
      MyOM->setVerbosity( Belos::Warnings|Belos::Debug);

      ghost_densemat_t* v = (ghost_densemat_t*)vec1_;
      Teuchos::RCP<phist::GhostMV> ivec = phist::rcp(v,false);

      // test the multivector and its adapter
      bool berr=Belos::TestMultiVecTraits<ST,phist::GhostMV>(MyOM,ivec);
      ASSERT_EQ(true,berr);
    }
  }
#endif
