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


TEST_F(CLASSNAME,fused_spmv__wdotw)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    std::vector<_MT_> v2dotv2(nvec_);
    std::vector<_MT_> v2nrm(nvec_);
    std::vector<_MT_> v3nrm(nvec_);

    // actually do y=Ax with ||y||
    SUBR(fused_spmv_mvdot)(alpha,A_,vec1_,beta,vec2_,&v2dotv2[0],NULL,&iflag_);
    ASSERT_EQ(0,iflag_);
    for (int i=0;i<nvec_;i++) v2nrm[i]=mt::sqrt(mt::real(v2dotv2[i]));

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check vnorm = ||y||
    SUBR(mvec_norm2)(vec3_,&v3nrm[0],&iflag_);
    for(int i = 0; i < nvec_; i++)
    {
      ASSERT_NEAR(v2nrm[i], v3nrm[i], 100*VTest::releps(vec3_));
    }
}

TEST_F(CLASSNAME,fused_spmv_vdotw)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    std::vector<_ST_> v2dotv1(nvec_);
    std::vector<_ST_> v3dotv1(nvec_);

    // actually do y=Ax with ||y||
    SUBR(fused_spmv_mvdot)(alpha,A_,vec1_,beta,vec2_,NULL,&v2dotv1[0],&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check vdot = x[i]' y[i]
    SUBR(mvec_dot_mvec)(vec3_,vec1_,&v3dotv1[0],&iflag_);
    for(int i = 0; i < nvec_; i++)
    {
      ASSERT_NEAR(st::real(v2dotv1[i]), st::real(v3dotv1[i]), 100*mt::sqrt(VTest::releps(vec3_)));
      ASSERT_NEAR(st::imag(v2dotv1[i]), st::imag(v3dotv1[i]), 100*mt::sqrt(VTest::releps(vec3_)));
    }
}

#if 1
#warning "these tests have to be updated to the new fused kernel interface!"
#else

TEST_F(CLASSNAME,fused_spmv_mvdot_mvadd)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    std::vector<_ST_> v2dotv1(nvec_);
    std::vector<_ST_> v3dotv1(nvec_);
    std::vector<_MT_> v2dotv2(nvec_);
    std::vector<_MT_> v3dotv3(nvec_);

    // actually do y=Ax and z=y-(AX)=0, and compute xdoty and ydoty
    SUBR(fused_spmv__mvdot_mvadd)(alpha,A_,vec1_,beta,vec2_,gamma,delta,vec3_,&v12dot[0],&v2nrm[0],&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check vdot = y[i]' x[i]
    SUBR(mvec_dot_mvec)(vec3_,vec1_,&v3dotv1[0],&iflag_);
    for(int i = 0; i < nvec_; i++)
    {
      ASSERT_NEAR(st::real(v2dotv1[i]), st::real(v3dotv1[i]), 100*mt::sqrt(VTest::releps(vec3_)));
      ASSERT_NEAR(st::imag(v2dotv1[i]), st::imag(v3dotv1[i]), 100*mt::sqrt(VTest::releps(vec3_)));
    }

    // check vnorm = ||y||
    SUBR(mvec_norm2)(vec3_,&v3nrm[0],&iflag_);
    for(int i = 0; i < _NV_; i++)
    {
      ASSERT_NEAR(v2nrm[i], v3nrm[i], 100*VTest::releps(vec3_));
    }
}



TEST_F(CLASSNAME,fused_spmv_wTw)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    // actually do y=Ax with yTy
    SUBR(fused_spmv_TROET__mvecT_times_mvec_self)(alpha,A_,vec1_,beta,vec2_,mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check yTy == y^T * y
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec3_,st::zero(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat2_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::sqrt(VTest::releps(vec3_)));
}


TEST_F(CLASSNAME,fused_spmv_wTv)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    // actually do y=Ax with yTx
    SUBR(fused_spmv_TROET__mvecT_times_mvec_other)(alpha,A_,vec1_,beta,vec2_,mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check yTx == y^T * x
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec1_,st::zero(),mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat2_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::sqrt(VTest::releps(vec3_)));
}


TEST_F(CLASSNAME,fused_spmv_TROET__mvTmv)
{
    int stride = 1;
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    // actually do y=Ax with yTx and yTy
    SUBR(fused_spmv_TROET__mvecT_times_mvec_both)(alpha,A_,vec1_,beta,vec2_,mat1_,mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check yTx == y^T * x
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec1_,st::zero(),mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat3_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::sqrt(VTest::releps(vec3_)));

    // check yTy == y^T * y
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec3_,st::zero(),mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat2_vp_,mat3_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::sqrt(VTest::releps(vec3_)));
}
#endif
