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


// test that we correctly get w'w and w=alpha*A*v+beta*w, and that it is ok to
// set the vdotw argument to NULL
TEST_F(CLASSNAME,fused_spmv_mvdot)
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
    
    std::vector<_ST_> dot_xy(nvec_);
    std::vector<_ST_> dot_yy(nvec_);
    std::vector<_ST_> dot_xy_ref(nvec_);
    std::vector<_ST_> dot_yy_ref(nvec_);
    
    // for checking y = alpha*A * x + beta*y
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(mvec_dot_mvec)(vec3_,vec1_,&dot_xy_ref[0],&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_dot_mvec)(vec3_,vec3_,&dot_yy_ref[0],&iflag_);
    ASSERT_EQ(0,iflag_);

    // actually do y=Ax with y'y and x'y
    SUBR(fused_spmv_mvdot)(alpha,A_,vec1_,beta,vec2_,&dot_yy[0],&dot_xy[0],&iflag_);    
    ASSERT_EQ(0,iflag_);

    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    for (int i=0;i<nvec_;i++)
    {
      ASSERT_NEAR(st::real(dot_xy[i]), st::real(dot_xy_ref[i]), 100*VTest::releps(vec3_));
      ASSERT_NEAR(st::imag(dot_xy[i]), st::imag(dot_xy_ref[i]), 100*VTest::releps(vec3_));
      ASSERT_NEAR(st::real(dot_yy[i]), st::real(dot_yy_ref[i]), 100*VTest::releps(vec3_));
      ASSERT_NEAR(st::imag(dot_yy[i]), st::imag(dot_yy_ref[i]), 100*VTest::releps(vec3_));
    }
}

TEST_F(CLASSNAME,fused_spmv_mvTmv)
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

    // actually do y=Ax with xTy and yTy
    SUBR(fused_spmv_mvTmv)(alpha,A_,vec1_,beta,vec2_,mat1_,mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check y = A * x
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride,vflag_), 100*VTest::releps(vec3_));

    // check xTy == x^T * y
    SUBR(mvecT_times_mvec)(st::one(),vec1_,vec3_,st::zero(),mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat2_vp_,mat3_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::sqrt(VTest::releps(vec3_)));

    // check yTy == y^T * y
    SUBR(mvecT_times_mvec)(st::one(),vec3_,vec3_,st::zero(),mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sdMat_from_device)(mat3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(mat1_vp_,mat3_vp_,nvec_,nvec_,m_lda_,stride,mflag_), 100*mt::sqrt(VTest::releps(vec3_)));
}

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
    
    // compute  v3=alpha*A*v1+beta*v2
    // and then v2=alpha*A*v1+beta*v2 'fused' v3=v3-v2(=0)
    
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(fused_spmv_mvdot_mvadd)(alpha,A_,vec1_,beta,vec2_,
        st::one(), -st::one(), vec3_, NULL,NULL,&iflag_);    
    ASSERT_EQ(0,iflag_);

    ASSERT_NEAR(mt::one(), MvecEqual(vec3_,st::zero()), 100*mt::eps());

}



