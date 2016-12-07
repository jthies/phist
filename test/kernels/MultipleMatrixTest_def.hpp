#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

using namespace phist::testing;

/*
#ifndef CONCAT
#define CONCAT(a,b) a ## b
#endif
#ifdef MTEST_NAME
#undef MTEST_NAME
#endif
#define VMTEST CONCAT(CLASSNAME,_VmTest) 

class VMTEST
: public KernelTestWithVectors<_ST_,_M_,_NV_,_USE_VIEWS_,3>
{
  VMTEST()
  {
  }
  ~VMTEST()
  {
  }
};
*/

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,3>
{

  public:
  
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,3> VTest;
  typedef TestWithType< _ST_ > ST_Test;
  typedef TestWithType< _MT_ > MT_Test;

  static void SetUpTestCase()
  {
    ST_Test::SetUpTestCase();
    KernelTest::SetUpTestCase();
    if (!typeImplemented_) return;
    int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
    // initialize all of the row functions that we use in this class
    ghost_gidx gdim[2];
    gdim[0]=_N_;
    gdim[1]=_N_;
    // initialize rowFuncs
    iflag_=phist::testing::PHIST_TG_PREFIX(hpd_tridiag)(-1,NULL,gdim,NULL,NULL);
    ASSERT_EQ(0,iflag_);
    iflag_=phist::testing::PHIST_TG_PREFIX(left_shift)(-1,NULL,gdim,NULL,NULL);
    ASSERT_EQ(0,iflag_);
    iflag_=phist::testing::PHIST_TG_PREFIX(left_shift_perio)(-1,NULL,gdim,NULL,NULL);
    ASSERT_EQ(0,iflag_);
    iflag_=phist::testing::PHIST_TG_PREFIX(right_shift)(-1,NULL,gdim,NULL,NULL);
    ASSERT_EQ(0,iflag_);
    iflag_=phist::testing::PHIST_TG_PREFIX(right_shift_perio)(-1,NULL,gdim,NULL,NULL);
    ASSERT_EQ(0,iflag_);
    iflag_=phist::testing::PHIST_TG_PREFIX(idfunc)(-1,NULL,gdim,NULL,NULL);
    ASSERT_EQ(0,iflag_);

    // create the main matrix, maps and context
    SUBR(sparseMat_create_fromRowFunc)(&A_,comm_,_N_,_N_,3,&phist::testing::PHIST_TG_PREFIX(hpd_tridiag),NULL,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_const_map_ptr range_map=NULL;
    SUBR(sparseMat_get_range_map)(A_,&range_map,&iflag_);
    ASSERT_EQ(0,iflag_);
        
    // initialize the base class with the context/maps
    KernelTestWithMap<_N_>::SetUpTestCaseWithMap(range_map);

    SUBR(sparseMat_get_context)(A_,&context_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // the other three matrices have a subset of the pattern of A each, so we should be able
    // to construct them with the same context:
    SUBR(sparseMat_create_fromRowFuncAndContext)(&I_,context_,1,&phist::testing::PHIST_TG_PREFIX(idfunc),NULL,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sparseMat_create_fromRowFuncAndContext)(&Al_,context_,1,&phist::testing::PHIST_TG_PREFIX(left_shift),NULL,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sparseMat_create_fromRowFuncAndContext)(&Ar_,context_,1,&phist::testing::PHIST_TG_PREFIX(right_shift),NULL,&iflag_);
    ASSERT_EQ(0,iflag_);

    // uninitialize rowFuncs
    iflag_=phist::testing::PHIST_TG_PREFIX(idfunc)(-2,NULL,NULL,NULL,NULL);
    iflag_=phist::testing::PHIST_TG_PREFIX(hpd_tridiag)(-2,NULL,NULL,NULL,NULL);
    iflag_=phist::testing::PHIST_TG_PREFIX(right_shift)(-2,NULL,NULL,NULL,NULL);
    iflag_=phist::testing::PHIST_TG_PREFIX(left_shift)(-2,NULL,NULL,NULL,NULL);

    VTest::SetUpTestCase();
  }

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    if (!typeImplemented_) return;
    VTest::SetUp();
    // initialize vec1=[[1:N]', [N+1:2*N]',...] and vec2=vec1
    int dim[2]; dim[0]=_N_; dim[1]=_NV_;
    SUBR(mvec_put_func)(vec1_,&PHIST_TG_PREFIX(mvec123func),dim,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    VTest::TearDown();
  }

  static void TearDownTestCase()
  {
    if (!typeImplemented_) return;
    VTest::TearDownTestCase();
    if (A_) SUBR(sparseMat_delete)(A_,&iflag_);
    if (I_) SUBR(sparseMat_delete)(I_,&iflag_);
    if (Al_) SUBR(sparseMat_delete)(Al_,&iflag_);
    if (Ar_) SUBR(sparseMat_delete)(Ar_,&iflag_);
    ASSERT_EQ(0,iflag_);
  }

protected:

  //! A is tridiagonal, I_ the identity matrix,
  //! Al and Au have only a unit  lower and upper
  //! diagonal, respectively.
  static TYPE(sparseMat_ptr) A_,I_,Al_,Ar_;
  
  //! the context of A_
  static phist_const_context_ptr context_;
  

};

TYPE(sparseMat_ptr) CLASSNAME :: A_;
TYPE(sparseMat_ptr) CLASSNAME :: Al_;
TYPE(sparseMat_ptr) CLASSNAME :: Ar_;
TYPE(sparseMat_ptr) CLASSNAME :: I_;
phist_const_context_ptr CLASSNAME::context_;

TEST_F(CLASSNAME,I_works)
{
  if (problemTooSmall_ || !typeImplemented_) return;
  // test X==I*X
  SUBR(sparseMat_times_mvec)(st::one(),I_,vec1_,st::zero(),vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);
  ASSERT_NEAR(mt::one(),MvecsEqual(vec1_,vec3_),mt::eps());
}

TEST_F(CLASSNAME,Al_works)
{
  if (problemTooSmall_ || !typeImplemented_) return;
  // test X==Al*X
  SUBR(sparseMat_times_mvec)(st::one(),Al_,vec1_,st::zero(),vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_from_device)(vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_from_device)(vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);
  // vec3(0:nloc-2,:) should be the same as vec2(1:nloc-1,:)
  ASSERT_NEAR(mt::one(),ArraysEqual(&(vec2_vp_[VIDX(1,0,lda_)]),&(vec3_vp_[VIDX(0,0,lda_)]),nloc_-1,_NV_,lda_,stride_,vflag_),mt::eps());
  //vec3(nloc-1,:)=vec2(0,:) on the next process, or 0 if rank==np-1
  for (int i=0; i<_NV_; i++)
  {
    _ST_ v2=vec2_vp_[VIDX(nloc_-1,0,lda_)];
    _ST_ v3=vec3_vp_[VIDX(nloc_-1,0,lda_)];
    _ST_ v_expect=((int)st::real(v2)%(_N_))? v2+st::one(): st::zero();
    ASSERT_REAL_EQ(st::real(v3),st::real(v_expect));
    ASSERT_REAL_EQ(st::imag(v3),st::imag(v_expect));
  }
}

TEST_F(CLASSNAME,Ar_works)
{
  if (problemTooSmall_ || !typeImplemented_) return;
  // test X==Ar*X
  SUBR(sparseMat_times_mvec)(st::one(),Ar_,vec1_,st::zero(),vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_from_device)(vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_from_device)(vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);
  // vec2(0:nloc-1,:) should be the same as vec3(1:nloc,:)
  ASSERT_NEAR(mt::one(),ArraysEqual(&(vec3_vp_[VIDX(1,0,lda_)]),&(vec2_vp_[VIDX(0,0,lda_)]),nloc_-1,_NV_,lda_,stride_,vflag_),mt::eps());
  //vec3(0,:)=vec2(nloc-1,:) on the previous process, or 0 if rank==0
  for (int i=0; i<_NV_; i++)
  {
    _ST_ v2=vec2_vp_[VIDX(0,0,lda_)];
    _ST_ v3=vec3_vp_[VIDX(0,0,lda_)];
    _ST_ v_expect=v2-st::one();
    ASSERT_REAL_EQ(st::real(v3),st::real(v_expect));
    ASSERT_REAL_EQ(st::imag(v3),st::imag(v_expect));
  }
}

//this test doesn't work in complex arithmetic at the moment
// because of assumptions on the A matrix, we disable it for 
// the moment
#ifdef IS_COMPLEX
TEST_F(CLASSNAME,DISABLED_A_works)
#else
TEST_F(CLASSNAME,A_works)
#endif
{
  if (problemTooSmall_ || !typeImplemented_) return;
  // test X==Ar*X
  SUBR(sparseMat_times_mvec)(st::one(),A_,vec1_,st::zero(),vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_from_device)(vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_from_device)(vec3_,&iflag_);
  ASSERT_EQ(0,iflag_);

  // this matrix is like a second derivative on the inner points, so it gives 0 unless i=0 or _N_-1.
  // For those we explicitly subtract the expected value, then compare everything to 0 
  if (mpi_rank_==0)
  {
    for (int i=0; i<_NV_; i++)
    {
      vec3_vp_[VIDX(0,i,lda_)] -= vec2_vp_[VIDX(0,i,lda_)] - (_ST_)0.5*vec2_vp_[VIDX(1,i,lda_)];
    }
  }
  if (mpi_rank_==mpi_size_-1)
  {
    for (int i=0; i<_NV_; i++)
    {
      vec3_vp_[VIDX(nloc_-1,i,lda_)] -= vec2_vp_[VIDX(nloc_-1,i,lda_)] - (_ST_)0.5*vec2_vp_[VIDX(nloc_-2,i,lda_)];
    }
  }
  ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,_NV_,lda_,stride_,st::zero(),vflag_));
}
