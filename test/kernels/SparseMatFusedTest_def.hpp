#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,3>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_,_USE_VIEWS_>
{

public:
  
  typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,_USE_VIEWS_,3> VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_,_USE_VIEWS_> MTest;
  typedef TestWithType< _MT_ > MT_Test;

  static void SetUpTestCase()
  {
    SparseMatTest::SetUpTestCase();
    VTest::SetUpTestCase();
    MTest::SetUpTestCase();
  }

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    SparseMatTest::SetUp();
    VTest::SetUp();
    MTest::SetUp();
    
    haveMat_ = (A_ != NULL);
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    MTest::TearDown();
    VTest::TearDown();
    SparseMatTest::TearDown();
  }

  static void TearDownTestCase()
  {
    MTest::TearDownTestCase();
    VTest::TearDownTestCase();
    SparseMatTest::TearDownTestCase();
  }

protected:
  bool haveMat_;
};


TEST_F(CLASSNAME,sparseMat_times_mvec_fused_norm2)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::rand();
    _ST_ beta = st::rand();
    std::vector<_MT_> v2nrm(nvec_);
    std::vector<_MT_> v3nrm(nvec_);

    // actually do y=Ax with ||y||
    SUBR(sparseMat_times_mvec_fused_norm2)(alpha,A_,vec1_,beta,vec2_,&v2nrm[0],&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*mt::eps());

    // check vnorm = ||y||
    SUBR(mvec_norm2)(vec3_,&v3nrm[0],&iflag_);
    for(int i = 0; i < nvec_; i++)
    {
      ASSERT_NEAR(v2nrm[i], v3nrm[i], 100*mt::eps());
    }
}

TEST_F(CLASSNAME,sparseMat_times_mvec_fused_dot)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::rand();
    _ST_ beta = st::rand();
    std::vector<_ST_> v12dot(nvec_);
    std::vector<_ST_> v13dot(nvec_);

    // actually do y=Ax with ||y||
    SUBR(sparseMat_times_mvec_fused_dot)(alpha,A_,vec1_,beta,vec2_,&v12dot[0],&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*mt::eps());

    // check vdot = y[i]' x[i]
    SUBR(mvec_dot_mvec)(vec3_,vec1_,&v13dot[0],&iflag_);
    for(int i = 0; i < nvec_; i++)
    {
      ASSERT_NEAR(st::real(v12dot[i]), st::real(v13dot[i]), 100*mt::eps());
      ASSERT_NEAR(st::imag(v12dot[i]), st::imag(v13dot[i]), 100*mt::eps());
    }
}

TEST_F(CLASSNAME,sparseMat_times_mvec_fused_dot_norm2)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::rand();
    _ST_ beta = st::rand();
    std::vector<_ST_> v12dot(nvec_);
    std::vector<_ST_> v13dot(nvec_);
    std::vector<_MT_> v2nrm(nvec_);
    std::vector<_MT_> v3nrm(nvec_);

    // actually do y=Ax with ||y||
    SUBR(sparseMat_times_mvec_fused_dot_norm2)(alpha,A_,vec1_,beta,vec2_,&v12dot[0],&v2nrm[0],&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*mt::eps());

    // check vdot = y[i]' x[i]
    SUBR(mvec_dot_mvec)(vec3_,vec1_,&v13dot[0],&iflag_);
    for(int i = 0; i < nvec_; i++)
    {
      ASSERT_NEAR(st::real(v12dot[i]), st::real(v13dot[i]), 100*mt::eps());
      ASSERT_NEAR(st::imag(v12dot[i]), st::imag(v13dot[i]), 100*mt::eps());
    }

    // check vnorm = ||y||
    SUBR(mvec_norm2)(vec3_,&v3nrm[0],&iflag_);
    for(int i = 0; i < _NV_; i++)
    {
      ASSERT_NEAR(v2nrm[i], v3nrm[i], 100*mt::eps());
    }
}



TEST_F(CLASSNAME,sparseMat_times_mvec_fused_mvecT_times_mvec_self)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::rand();
    _ST_ beta = st::rand();

    // actually do y=Ax with yTy
    SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_self)(alpha,A_,vec1_,beta,vec2_,mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*mt::eps());

    // check yTy == y^T * y
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec3_,st::zero(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat2_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::eps());
}


TEST_F(CLASSNAME,sparseMat_times_mvec_fused_mvecT_times_mvec_other)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::rand();
    _ST_ beta = st::rand();

    // actually do y=Ax with yTx
    SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_other)(alpha,A_,vec1_,beta,vec2_,mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*mt::eps());

    // check yTx == y^T * x
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec1_,st::zero(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat2_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::eps());
}


TEST_F(CLASSNAME,sparseMat_times_mvec_fused_mvecT_times_mvec_both)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::rand();
    _ST_ beta = st::rand();

    // actually do y=Ax with yTx and yTy
    SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_both)(alpha,A_,vec1_,beta,vec2_,mat1_,mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*mt::eps());

    // check yTx == y^T * x
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec1_,st::zero(),mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat3_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::eps());

    // check yTy == y^T * y
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec3_,st::zero(),mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat2_vp_,mat3_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::eps());
}

