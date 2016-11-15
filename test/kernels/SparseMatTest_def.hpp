#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

#if MATNAME==MATNAME_BENCH3D_8_A1
# if defined(IS_COMPLEX)||!defined(IS_DOUBLE)
#   define DONT_INSTANTIATE
# endif
#endif


#ifndef DONT_INSTANTIATE

using namespace phist::testing;


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
    int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
    SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
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


  void test_sparseMat_times_mvec_vadd_mvec(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ shifts[_NV_], _ST_ beta)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts, vec1_, beta, vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    _ST_ alpha_shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      alpha_shifts[i] = alpha*shifts[i];

    SUBR(sparseMat_times_mvec)(alpha, A, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_vadd_mvec)(alpha_shifts, vec1_, st::one(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    ASSERT_NEAR(mt::one(), MvecsEqual(vec2_,vec3_,mt::one()),1000*mt::eps());
  }


  void test_sparseMat_times_mvec_communicate(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ beta)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_communicate)(A, vec1_, &iflag_);
    iflag_ = PHIST_SPMVM_ONLY_LOCAL;
    SUBR(sparseMat_times_mvec)(alpha, A, vec1_, beta, vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec)(alpha, A, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    ASSERT_NEAR(mt::one(), MvecsEqual(vec2_,vec3_), sqrt(mt::eps()));
  }


  void test_sparseMat_times_mvec_vadd_mvec_communicate(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ shifts[_NV_], _ST_ beta)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_communicate)(A, vec1_, &iflag_);
    iflag_ = PHIST_SPMVM_ONLY_LOCAL;
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts, vec1_, beta, vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    ASSERT_NEAR(mt::one(), MvecsEqual(vec2_,vec3_), sqrt(mt::eps()));
  }


#if _NV_ > 1
  void test_sparseMat_times_mvec_with_views(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ beta, int imin, int imax)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    // create views
    TYPE(mvec_ptr) vec1_view = NULL;
    SUBR(mvec_view_block)(vec1_, &vec1_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);
    TYPE(mvec_ptr) vec2_view = NULL;
    SUBR(mvec_view_block)(vec2_, &vec2_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec)(alpha, A, vec1_view, beta, vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);

    // make sure nothing changed outside of viewed block
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0, iflag_);
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,imin,lda_,stride_,vflag_));
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imax+1,lda_),vec3_vp_+VIDX(0,imax+1,lda_),nloc_,nvec_-imax-1,lda_,stride_,vflag_));

    // calculation for full block as reference
    SUBR(sparseMat_times_mvec)(alpha, A, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0, iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imin,lda_),vec3_vp_+VIDX(0,imin,lda_),nloc_,imax-imin+1,lda_,stride_,vflag_), mt::sqrt(VTest::releps()));

    // delete view
    SUBR(mvec_delete)(vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_delete)(vec1_view, &iflag_);
    ASSERT_EQ(0, iflag_);
  }

  void test_sparseMat_times_mvec_vadd_mvec_with_views(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ shifts[_NV_], _ST_ beta, int imin, int imax)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    // create a view
    TYPE(mvec_ptr) vec1_view = NULL;
    SUBR(mvec_view_block)(vec1_, &vec1_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);
    TYPE(mvec_ptr) vec2_view = NULL;
    SUBR(mvec_view_block)(vec2_, &vec2_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts+imin, vec1_view, beta, vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0, iflag_);
    // make sure nothing changed outside of viewed block
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,imin,lda_,stride_,vflag_));
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imax+1,lda_),vec3_vp_+VIDX(0,imax+1,lda_),nloc_,nvec_-imax-1,lda_,stride_,vflag_));

    // calculation for full block as reference
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0, iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imin,lda_),vec3_vp_+VIDX(0,imin,lda_),nloc_,imax-imin+1,lda_,stride_,vflag_), mt::sqrt(VTest::releps(vec3_)));

    // delete view
    SUBR(mvec_delete)(vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_delete)(vec1_view, &iflag_);
    ASSERT_EQ(0, iflag_);
  }
#endif
  

protected:

  _MT_ const_row_sum_test(TYPE(sparseMat_ptr) A)
  {
    if (typeImplemented_ && !problemTooSmall_ && haveMat_)
    {
      _ST_ val = st::prand();
      global_sum(&val,1,mpi_comm_);
      SUBR(mvec_put_value)(vec1_,val,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A,vec1_,st::zero(),vec2_,&iflag_);
      if (iflag_) return (_MT_)iflag_;

      PrintVector(PHIST_DEBUG,"input to spMVM",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"result of spMVM",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      return MvecEqual(vec2_,val);
    }
    return mt::one();
  }

  bool haveMat_;
};

  TEST_F(CLASSNAME, read_matrices) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_TRUE(AssertNotNull(A_));
    
      // test that the global number of rows/cols is correct in the objects
      phist_gidx gnrows, gncols;
      SUBR(sparseMat_global_nrows)(A_,&gnrows,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(gnrows,nglob_);
      SUBR(sparseMat_global_ncols)(A_,&gncols,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_EQ(gncols,nglob_);
    }
  }

#if MATNAME == MATNAME_spzero
  TEST_F(CLASSNAME, A0_times_mvec) 
  {
    if (typeImplemented_ && !problemTooSmall_ && haveMat_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec2_,0.0));
    }
  }
#endif

#if MATNAME == MATNAME_speye
  TEST_F(CLASSNAME, A1_times_mvec)
  {
    if (typeImplemented_ && !problemTooSmall_ && haveMat_)
    {
      ST alpha, beta;
      //I*X=X?
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec1_,vec2_));

      //alpha*I*X=alpha*X?
      alpha = st::prand();
      beta=st::zero();
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_scale)(vec1_,alpha,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec1_,vec2_));

      //0*I*X+beta*Y = beta*Y? 
      alpha=st::zero(); 
      beta=st::prand();
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
#endif
      SUBR(mvec_random)(vec1_,&iflag_); 
      SUBR(mvec_random)(vec2_,&iflag_); 

      PrintVector(PHIST_DEBUG,"input to spMVM",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"result of spMVM",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);

      // v3=beta*v2 
      SUBR(mvec_add_mvec)(beta,vec2_,st::zero(),vec3_,&iflag_); 
      ASSERT_EQ(0,iflag_); 
      // v2 = 0*v1 + beta*v2 (=v3) 
      SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_); 
      ASSERT_EQ(0,iflag_); 

      PrintVector(PHIST_DEBUG,"expected",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"result of spMVM with beta!=0",vec3_vp_,nloc_,lda_,stride_,mpi_comm_);

      ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec2_,vec3_));

      //I*X+beta*Y = X+beta*Y?
      alpha = st::one();
      beta = st::prand();
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
#endif
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);

      PrintVector(PHIST_DEBUG,"input",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"output, before",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);

      //v3=v1+beta*v2
      SUBR(mvec_to_mvec)(vec1_,vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_scale)(vec3_, alpha, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // v2 = v1 + beta*v2 (=alpha*v1+v3)
      SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      PrintVector(PHIST_DEBUG,"result with beta!=0",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"expected",vec3_vp_,nloc_,lda_,stride_,mpi_comm_);
\
      ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_,mt::one()),1000*mt::eps());

      //alpha*I*X+beta*Y = alpha*X+beta*Y?
      alpha = st::prand();
      beta = st::prand();
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
#endif
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);

      PrintVector(PHIST_DEBUG,"input",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"output, before",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);

       // v3=alpha*v1+beta*v2
      SUBR(mvec_to_mvec)(vec1_,vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_scale)(vec3_, alpha, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // v2 = alpha*v1 + beta*v2 (=alpha*v1+v3)
      SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      PrintVector(PHIST_DEBUG,"result with beta!=0",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      PrintVector(PHIST_DEBUG,"expected",vec3_vp_,nloc_,lda_,stride_,mpi_comm_);

      ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_,mt::one()),1000*mt::eps());
    }
  }

#if(_NV_>1)
  TEST_F(CLASSNAME, A1_times_mvec_using_two_views_of_the_same_vec)
  {
    if (typeImplemented_ && !problemTooSmall_ && A_!=NULL)
    {
      ST alpha, beta;
      //I*X=X?
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      
      TYPE(mvec_ptr) v_in=NULL, v_out=NULL;
      phist_lidx nv =  _NV_/2;
      phist_lidx offs = _NV_%2;
      SUBR(mvec_view_block)(vec1_,&v_in,offs,offs+nv-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_view_block)(vec1_,&v_out,offs+nv,offs+2*nv-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      
     // _ST_ *v_in_vp, *v_out_vp;

// we could extract the pointers from the view, but this way we
// have more control and make sure the data is where it should be:
     //v_in_vp = vec1_vp_+VIDX(0,offs,lda_);
     //v_out_vp = vec1_vp_+VIDX(0,offs+nv,lda_);
      PrintVector(PHIST_DEBUG,"random",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);

      SUBR(sparseMat_times_mvec)(st::one(),A_,v_in,st::zero(),v_out,&iflag_);
      ASSERT_EQ(0,iflag_);

      PHIST_SOUT(PHIST_DEBUG,"with two identical blocks [%d..%d] and [%d..%d]:\n",
        offs, offs+nv-1,offs+nv,offs+2*nv-1);
      PrintVector(PHIST_DEBUG,"with two identical blocks",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);

      ASSERT_REAL_EQ(mt::one(),MvecsEqual(v_in,v_out));

      SUBR(mvec_delete)(v_out,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(v_in,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    else
    {
      PHIST_SOUT(PHIST_INFO,"test skipped\n");
    }
  }
#endif

  TEST_F(CLASSNAME,A1_fromRowFunc)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    TYPE(sparseMat_ptr) A=NULL;
    iflag_=PHIST_SPARSEMAT_QUIET;
    SUBR(sparseMat_create_fromRowFunc)(&A,comm_,_N_,_N_,1,&PHIST_TG_PREFIX(idfunc),NULL,&iflag_);
    // use EXPECT from now on so that the matrix gets deleted even if something fails
    ASSERT_EQ(0,iflag_);

    // the domain map of this matrix should be the defaultMap_, check this and construct vectors for the spMVM
    phist_const_map_ptr domainA=NULL;
    SUBR(sparseMat_get_domain_map)(A,&domainA,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_maps_compatible(defaultMap_,domainA,&iflag_);
    EXPECT_EQ(0,iflag_);
    // and it has to be compatible with the (domain) map of A_ that defines the vectors vec1/2_
    phist_maps_compatible(map_,domainA,&iflag_);
    EXPECT_TRUE(iflag_>=0);
    bool other_map=(iflag_>0);
    TYPE(mvec_ptr) vec1=vec1_, vec2=vec2_;
    if (other_map)
    {
      PHISTTEST_MVEC_CREATE(&vec1,domainA,nvec_,&iflag_);
      EXPECT_EQ(0,iflag_);
      PHISTTEST_MVEC_CREATE(&vec2,domainA,nvec_,&iflag_);
      EXPECT_EQ(0,iflag_);
    }
    
    // check that AX=X
    SUBR(mvec_random)(vec1,&iflag_);
    EXPECT_EQ(0,iflag_);
    SUBR(mvec_put_value)(vec2,(_ST_)99.9,&iflag_);
    EXPECT_EQ(0,iflag_);
    SUBR(sparseMat_times_mvec)(st::one(),A,vec1,st::zero(),vec2,&iflag_);
    EXPECT_EQ(0,iflag_);
    EXPECT_NEAR(1.0,MvecsEqual(vec1,vec2),100*mt::eps());
    if (other_map)
    {
      SUBR(mvec_delete)(vec1,&iflag_);
      EXPECT_EQ(0,iflag_);
      SUBR(mvec_delete)(vec2,&iflag_);
      EXPECT_EQ(0,iflag_);
    }
    SUBR(sparseMat_delete)(A,&iflag_);
    EXPECT_EQ(0,iflag_);
  }

  TEST_F(CLASSNAME,A1_fromRowFuncAndContext)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    TYPE(sparseMat_ptr) A=NULL;
    iflag_=PHIST_SPARSEMAT_QUIET;
    SUBR(sparseMat_create_fromRowFuncAndContext)(&A,context_,1,&PHIST_TG_PREFIX(idfunc),NULL,&iflag_);
    ASSERT_EQ(0,iflag_);

    // check that AX=X
    SUBR(mvec_random)(vec1_,&iflag_);
    EXPECT_EQ(0,iflag_);
    SUBR(mvec_put_value)(vec2_,(_ST_)99.9,&iflag_);
    EXPECT_EQ(0,iflag_);
    SUBR(sparseMat_times_mvec)(st::one(),A,vec1_,st::zero(),vec2_,&iflag_);
    EXPECT_EQ(0,iflag_);
    EXPECT_NEAR(1.0,MvecsEqual(vec1_,vec2_),100*mt::eps());
    SUBR(sparseMat_delete)(A,&iflag_);
    EXPECT_EQ(0,iflag_);
  }

#endif // MATNAME_speye

#if MATNAME == MATNAME_sprandn
  TEST_F(CLASSNAME, A2_times_mvec)
    {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // we allow a tolerance here because the matrices may have errors in the
    // last digit and we can't get the test to pass otherwise.
    ASSERT_NEAR(mt::one(),const_row_sum_test(A_),100*VTest::releps());
    }
#endif

#if MATNAME == MATNAME_sprandn_nodiag
  TEST_F(CLASSNAME, A3_times_mvec)
    {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // we allow a tolerance here because the matrices may have errors in the
    // last digit and we can't get the test to pass otherwise.
    ASSERT_NEAR(mt::one(),const_row_sum_test(A_),100*VTest::releps());
    }
#endif

#if MATNAME == MATNAME_spshift
  TEST_F(CLASSNAME, shift_mvec)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      _ST_ alpha = st::one();
      _ST_ beta = st::zero();

      // create new map (which has correct order!)
      phist_map_ptr map;
      phist_map_create(&map, comm_, nglob_, &iflag_);
      ASSERT_EQ(0,iflag_);
      phist_lidx nloc = 0;
      phist_map_get_local_length(map,&nloc,&iflag_);
      ASSERT_EQ(0,iflag_);
      phist_gidx ilower = 0;
      phist_map_get_ilower(map,&ilower,&iflag_);
      ASSERT_EQ(0,iflag_);

      // create a vector with this map
      TYPE(mvec_ptr) orderedVec = NULL;
      PHISTTEST_MVEC_CREATE(&orderedVec, map, nvec_, &iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ *orderedVec_vp = NULL;
      phist_lidx lda = 0;
      SUBR(mvec_extract_view)(orderedVec, &orderedVec_vp, &lda, &iflag_);
      ASSERT_EQ(0,iflag_);

      // setup recognizable input
      for(int i = 0; i < nloc; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          orderedVec_vp[VIDX(i,j,lda)] = (_ST_)(ilower+i + j*nglob_);
        }
      }
      // copy to vec1_
      SUBR(mvec_to_mvec)(orderedVec, vec1_, &iflag_);
      ASSERT_EQ(0,iflag_);


      // apply our shift matrix
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A='shift', alpha=1, beta=0"<<std::endl;
#endif

      PrintVector(PHIST_DEBUG,"input",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);

      SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      PrintVector(PHIST_DEBUG,"output",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);

      // copy to orderedVec
      SUBR(mvec_to_mvec)(vec2_, orderedVec, &iflag_);
      ASSERT_EQ(0,iflag_);


      // check result
      for(int i = 0; i < nglob_; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          if( i >= ilower && i < ilower+nloc )
          {
            ASSERT_REAL_EQ((i+1)%nglob_ + j*nglob_,st::real(orderedVec_vp[VIDX(i-ilower,j,lda)]));
          }
          else
          {
            // also assert elements of other processes, s.t. all processes do the same number of asserts!
            ASSERT_REAL_EQ((i+1)%nglob_ + j*nglob_, (i+1)%nglob_ + j*nglob_);
          }
        }
      }

      SUBR(mvec_delete)(orderedVec,&iflag_);
      ASSERT_EQ(0,iflag_);
      phist_map_delete(map,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
#endif


#if MATNAME == MATNAME_sprandn
  TEST_F(CLASSNAME, A2_precalc_result)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      _ST_ alpha = st::one();
      _ST_ beta = st::zero();

      // create new map (which has correct order!)
      // The test below relies on the function mvec_to_mvec
      // to translate the vectors back and forth between
      // this map and the matrix' domain map. This is the
      // case with GHOST if the matrix is not repartitioned
      // right now (that is, the difference between the map
      // is a local permutation induced e.g. by sigma>1 in
      // the SELL-C-sigma format).
      phist_map_ptr map;
      phist_map_create(&map, comm_, nglob_, &iflag_);
      ASSERT_EQ(0,iflag_);
      phist_lidx nloc = 0;
      phist_map_get_local_length(map,&nloc,&iflag_);
      ASSERT_EQ(0,iflag_);
      phist_gidx ilower = 0;
      phist_map_get_ilower(map,&ilower,&iflag_);
      ASSERT_EQ(0,iflag_);

      // create a vector with this map
      TYPE(mvec_ptr) orderedVec = NULL;
      PHISTTEST_MVEC_CREATE(&orderedVec, map, nvec_, &iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ *orderedVec_vp = NULL;
      phist_lidx lda = 0;
      SUBR(mvec_extract_view)(orderedVec, &orderedVec_vp, &lda, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(mvec_from_device)(orderedVec, &iflag_);

      // setup recognizable input
      for(int i = 0; i < nloc; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          orderedVec_vp[VIDX(i,j,lda)] = (_ST_)(ilower+i + j*nglob_);
        }
      }
      
      SUBR(mvec_to_device)(orderedVec, &iflag_);
      
      // copy to vec1_
      SUBR(mvec_to_mvec)(orderedVec, vec1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      

      // apply our shift matrix
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A='rand', alpha=1, beta=0"<<std::endl;
#endif

      PrintVector(PHIST_DEBUG,"input",vec1_vp_,nloc_,lda_,stride_,mpi_comm_);

      SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);      
      ASSERT_EQ(0,iflag_);

      PrintVector(PHIST_DEBUG,"output",vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
 
#if _N_ == 25 && _NV_ == 1
#ifdef IS_COMPLEX
      _ST_ precalc_result[_N_*_NV_] = {
        _ST_(1.072230172746e+01,-1.031665649804e+00), 
        _ST_(2.131001265990e+01,-3.282812619721e+00), 
        _ST_(2.015641860460e+01,-1.541996705535e+01), 
        _ST_(2.171736316821e+01,6.518151979141e+00), 
        _ST_(2.615219127525e+00,-3.620215930236e+00), 
        _ST_(1.228901252787e+01,-4.528370474365e+00), 
        _ST_(-6.543281193088e+00,4.169322949310e+00), 
        _ST_(1.513919626135e+01,-3.968740334506e+00), 
        _ST_(3.010875585729e+01,-7.275935201946e+00), 
        _ST_(1.856357361375e+01,-6.082795919035e-01), 
        _ST_(4.884622853589e+00,-1.246056153190e+01), 
        _ST_(9.724002250734e+00,1.854295108448e+00), 
        _ST_(1.295234521942e+01,4.797517156124e+00), 
        _ST_(1.922176454031e+01,-6.554985875651e+00), 
        _ST_(9.937028268565e+00,-3.048436869856e+00), 
        _ST_(8.417823023929e+00,1.504269262034e+00), 
        _ST_(2.827516113199e+00,-6.272713645344e+00), 
        _ST_(1.542580758006e+01,-5.650405298966e+00), 
        _ST_(1.507749351511e+01,-4.585936779631e+00), 
        _ST_(1.644503601436e+01,-6.053982924770e+00), 
        _ST_(3.407669499649e+01,1.638938366859e+01), 
        _ST_(2.078494495946e+01,5.856000872126e-01), 
        _ST_(1.420728160837e+01,1.050663583365e+00), 
        _ST_(1.251387547773e+01,-4.875170138989e-01), 
        _ST_(3.325495196562e+00,-1.806246226644e+01) };
#else
      _ST_ precalc_result[_N_*_NV_] = {
        -0.12861779756316816,
        16.076787993876188,
        20.072300560349166,
        17.898973758721571,
        19.186951140368311,
        -52.388079187547561,
        12.117945192440743,
        11.417408852483819,
        18.080795442101312,
        13.972038117529888,
        4.3052769861597753,
        12.535799082066269,
        5.9149777593318920,
        8.3650228590087750E-002,
        10.824168487154179,
        44.098112302021931,
        80.769125851754040,
        10.243308300563234,
        -40.240733293758780,
        16.598578972374874,
        18.400975207367228,
        36.482850494569007,
        13.997556100937755,
        -8.1932230989204431,
        9.6675923086516082 };
#endif
#elif _N_ == 25 && _NV_ == 4
#ifdef IS_COMPLEX
      _ST_ precalc_result[_N_*_NV_] = {
        _ST_(1.072230172746e+01,-1.031665649804e+00), _ST_(3.572230172746e+01,-1.031665649804e+00), _ST_(6.072230172746e+01,-1.031665649804e+00), _ST_(8.572230172746e+01,-1.031665649804e+00),
        _ST_(2.131001265990e+01,-3.282812619721e+00), _ST_(4.631001265990e+01,-3.282812619721e+00), _ST_(7.131001265990e+01,-3.282812619721e+00), _ST_(9.631001265990e+01,-3.282812619721e+00),
        _ST_(2.015641860460e+01,-1.541996705535e+01), _ST_(4.515641860460e+01,-1.541996705535e+01), _ST_(7.015641860460e+01,-1.541996705535e+01), _ST_(9.515641860460e+01,-1.541996705535e+01),
        _ST_(2.171736316821e+01,6.518151979141e+00), _ST_(4.671736316821e+01,6.518151979141e+00), _ST_(7.171736316821e+01,6.518151979141e+00), _ST_(9.671736316821e+01,6.518151979141e+00),
        _ST_(2.615219127525e+00,-3.620215930236e+00), _ST_(2.761521912752e+01,-3.620215930236e+00), _ST_(5.261521912752e+01,-3.620215930236e+00), _ST_(7.761521912752e+01,-3.620215930236e+00),
        _ST_(1.228901252787e+01,-4.528370474365e+00), _ST_(3.728901252787e+01,-4.528370474365e+00), _ST_(6.228901252787e+01,-4.528370474365e+00), _ST_(8.728901252787e+01,-4.528370474365e+00),
        _ST_(-6.543281193088e+00,4.169322949310e+00), _ST_(1.845671880691e+01,4.169322949310e+00), _ST_(4.345671880691e+01,4.169322949310e+00), _ST_(6.845671880691e+01,4.169322949310e+00),
        _ST_(1.513919626135e+01,-3.968740334506e+00), _ST_(4.013919626135e+01,-3.968740334506e+00), _ST_(6.513919626135e+01,-3.968740334506e+00), _ST_(9.013919626135e+01,-3.968740334506e+00),
        _ST_(3.010875585729e+01,-7.275935201946e+00), _ST_(5.510875585729e+01,-7.275935201946e+00), _ST_(8.010875585729e+01,-7.275935201946e+00), _ST_(1.051087558573e+02,-7.275935201946e+00),
        _ST_(1.856357361375e+01,-6.082795919035e-01), _ST_(4.356357361375e+01,-6.082795919035e-01), _ST_(6.856357361375e+01,-6.082795919035e-01), _ST_(9.356357361375e+01,-6.082795919035e-01),
        _ST_(4.884622853589e+00,-1.246056153190e+01), _ST_(2.988462285359e+01,-1.246056153190e+01), _ST_(5.488462285359e+01,-1.246056153190e+01), _ST_(7.988462285359e+01,-1.246056153190e+01),
        _ST_(9.724002250734e+00,1.854295108448e+00), _ST_(3.472400225073e+01,1.854295108448e+00), _ST_(5.972400225073e+01,1.854295108448e+00), _ST_(8.472400225073e+01,1.854295108448e+00),
        _ST_(1.295234521942e+01,4.797517156124e+00), _ST_(3.795234521942e+01,4.797517156124e+00), _ST_(6.295234521942e+01,4.797517156124e+00), _ST_(8.795234521942e+01,4.797517156124e+00),
        _ST_(1.922176454031e+01,-6.554985875651e+00), _ST_(4.422176454031e+01,-6.554985875651e+00), _ST_(6.922176454031e+01,-6.554985875651e+00), _ST_(9.422176454031e+01,-6.554985875651e+00),
        _ST_(9.937028268565e+00,-3.048436869856e+00), _ST_(3.493702826857e+01,-3.048436869856e+00), _ST_(5.993702826857e+01,-3.048436869856e+00), _ST_(8.493702826857e+01,-3.048436869856e+00),
        _ST_(8.417823023929e+00,1.504269262034e+00), _ST_(3.341782302393e+01,1.504269262034e+00), _ST_(5.841782302393e+01,1.504269262034e+00), _ST_(8.341782302393e+01,1.504269262034e+00),
        _ST_(2.827516113199e+00,-6.272713645344e+00), _ST_(2.782751611320e+01,-6.272713645344e+00), _ST_(5.282751611320e+01,-6.272713645344e+00), _ST_(7.782751611320e+01,-6.272713645344e+00),
        _ST_(1.542580758006e+01,-5.650405298966e+00), _ST_(4.042580758006e+01,-5.650405298966e+00), _ST_(6.542580758006e+01,-5.650405298966e+00), _ST_(9.042580758006e+01,-5.650405298966e+00),
        _ST_(1.507749351511e+01,-4.585936779631e+00), _ST_(4.007749351511e+01,-4.585936779631e+00), _ST_(6.507749351511e+01,-4.585936779631e+00), _ST_(9.007749351511e+01,-4.585936779631e+00),
        _ST_(1.644503601436e+01,-6.053982924770e+00), _ST_(4.144503601436e+01,-6.053982924770e+00), _ST_(6.644503601436e+01,-6.053982924770e+00), _ST_(9.144503601436e+01,-6.053982924770e+00),
        _ST_(3.407669499649e+01,1.638938366859e+01), _ST_(5.907669499649e+01,1.638938366859e+01), _ST_(8.407669499649e+01,1.638938366859e+01), _ST_(1.090766949965e+02,1.638938366859e+01),
        _ST_(2.078494495946e+01,5.856000872126e-01), _ST_(4.578494495946e+01,5.856000872126e-01), _ST_(7.078494495946e+01,5.856000872126e-01), _ST_(9.578494495946e+01,5.856000872125e-01),
        _ST_(1.420728160837e+01,1.050663583365e+00), _ST_(3.920728160837e+01,1.050663583365e+00), _ST_(6.420728160837e+01,1.050663583365e+00), _ST_(8.920728160837e+01,1.050663583365e+00),
        _ST_(1.251387547773e+01,-4.875170138989e-01), _ST_(3.751387547773e+01,-4.875170138989e-01), _ST_(6.251387547773e+01,-4.875170138989e-01), _ST_(8.751387547773e+01,-4.875170138989e-01),
        _ST_(3.325495196562e+00,-1.806246226644e+01), _ST_(2.832549519656e+01,-1.806246226644e+01), _ST_(5.332549519656e+01,-1.806246226644e+01), _ST_(7.832549519656e+01,-1.806246226644e+01) };
#else
      _ST_ precalc_result[_N_*_NV_] = {
        -0.12861779756316816, 24.871382202436834, 49.871382202436841, 74.871382202436862,
        16.076787993876188, 41.076787993876181, 66.076787993876167, 91.076787993876181,
        20.072300560349166, 45.072300560349177, 70.072300560349163, 95.072300560349149,
        17.898973758721571, 42.898973758721567, 67.898973758721596, 92.898973758721581,
        19.186951140368311, 44.186951140368315, 69.186951140368308, 94.186951140368308,
        -52.388079187547561, -27.388079187547561, -2.3880791875475609, 22.611920812452411,
        12.117945192440743, 37.117945192440750, 62.117945192440757, 87.117945192440772,
        11.417408852483819, 36.417408852483824, 61.417408852483824, 86.417408852483831,
        18.080795442101312, 43.080795442101319, 68.080795442101319, 93.080795442101319,
        13.972038117529888, 38.972038117529891, 63.972038117529905, 88.972038117529905,
        4.3052769861597753, 29.305276986159761, 54.305276986159726, 79.305276986159711,
        12.535799082066269, 37.535799082066255, 62.535799082066234, 87.535799082066234,
        5.9149777593318920, 30.914977759331908, 55.914977759331904, 80.914977759331919,
        8.3650228590087750E-002, 25.083650228590084, 50.083650228590088, 75.083650228590102,
        10.824168487154179, 35.824168487154182, 60.824168487154182, 85.824168487154196,
        44.098112302021931, 69.098112302021960, 94.098112302021960, 119.09811230202197,
        80.769125851754040, 105.76912585175401, 130.76912585175404, 155.76912585175393,
        10.243308300563234, 35.243308300563228, 60.243308300563228, 85.243308300563228,
        -40.240733293758780, -15.240733293758730, 9.7592667062412914, 34.759266706241320,
        16.598578972374874, 41.598578972374867, 66.598578972374881, 91.598578972374881,
        18.400975207367228, 43.400975207367232, 68.400975207367239, 93.400975207367239,
        36.482850494569007, 61.482850494569007, 86.482850494569036, 111.48285049456905,
        13.997556100937755, 38.997556100937771, 63.997556100937771, 88.997556100937743,
        -8.1932230989204431, 16.806776901079530, 41.806776901079516, 66.806776901079530,
        9.6675923086516082, 34.667592308651606, 59.667592308651606, 84.667592308651621 };
#endif
#else
      PHIST_SOUT(PHIST_INFO,"skipping test, no data available\n");
      ST* precalc_result=NULL;
      SUBR(mvec_delete)(orderedVec,&iflag_);
      ASSERT_EQ(0,iflag_);
      phist_map_delete(map,&iflag_);
      ASSERT_EQ(0,iflag_);
      return;
#endif
      // copy to orderedVec
      SUBR(mvec_to_mvec)(vec2_, orderedVec, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      // download result
      SUBR(mvec_from_device)(orderedVec, &iflag_);

      // check result
      _MT_ err = 0;
      for(int i = 0; i < nloc; i++)
        for(int j = 0; j < nvec_; j++)
          err = std::max(err, st::abs( precalc_result[(ilower+i)*nvec_+j] - orderedVec_vp[VIDX(i,j,lda)] ));
      ASSERT_NEAR(err, mt::zero(), mt::sqrt(mt::eps()));

      SUBR(mvec_delete)(orderedVec,&iflag_);
      ASSERT_EQ(0,iflag_);
      phist_map_delete(map,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
#endif

#if MATNAME == MATNAME_sprandn
  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_only_scale)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::zero();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A_, shifts, beta);
  }
#endif

#if MATNAME == MATNAME_spzero
  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_only_vadd)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::one();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A_, shifts, beta);
  }
#endif

#if MATNAME == MATNAME_sprandn
  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_only_spmvm)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::zero();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A_, shifts, beta);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A_, shifts, beta);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_communicate_random)
  {
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_communicate(alpha, A_, beta);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_communicate_random)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_communicate(alpha, A_, shifts, beta);
  }
#endif


#if MATNAME == MATNAME_sprandn
#if _NV_ > 1
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_0_0)
  {
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A_, beta, 0, 0);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_0_0)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A_, shifts, beta, 0, 0);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_1_1)
  {
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A_, beta, 1, 1);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_1_1)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A_, shifts, beta, 1, 1);
  }
#endif


#if _NV_ > 3
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_1_2)
  {
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A_, beta, 1, 2);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_1_2)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A_, shifts, beta, 1, 2);
  }
#endif


#if _NV_ > 3
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_2_3)
  {
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A_, beta, 2, 3);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_2_3)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A_, shifts, beta, 2, 3);
  }
#endif


#if _NV_ > 4
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_1_4)
  {
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A_, beta, 1, 4);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_1_4)
  {
    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A_, shifts, beta, 1, 4);
  }
#endif


#if _NV_ >= 4
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_same_vec_views_0_1__2_3)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // random data
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0,iflag_);

    // save vec1_
    SUBR(mvec_add_mvec)(st::one(), vec1_, st::zero(), vec2_, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_to_mvec)(vec1_,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    // this is a sanity check to make sure the above give the same (mvec_add_mvec and mvec_to_mvec)
    ASSERT_REAL_EQ(1.0,MvecsEqual(vec1_,vec2_));
    ASSERT_REAL_EQ(1.0,MvecsEqual(vec2_,vec3_));

    // create views
    TYPE(mvec_ptr) vin = NULL;
    SUBR(mvec_view_block)(vec1_, &vin, 2, 3, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vout = NULL;
    SUBR(mvec_view_block)(vec1_, &vout, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vref = NULL;
    SUBR(mvec_view_block)(vec2_, &vref, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);

    // first generate reference data (safe calculation)
    SUBR(sparseMat_times_mvec)(alpha, A_, vin, beta, vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    // calculation (unsafe aliasing!)
    SUBR(sparseMat_times_mvec)(alpha, A_, vin, beta, vout, &iflag_);
    ASSERT_EQ(0,iflag_);

    // check vin
    SUBR(mvec_from_device)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_+VIDX(0,2,lda_),vec3_vp_+VIDX(0,2,lda_),nloc_,2,lda_,stride_,vflag_), mt::eps());
    // check result
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_,vec2_vp_,nloc_,2,lda_,stride_,vflag_), sqrt(mt::eps()));

    // delete views
    SUBR(mvec_delete)(vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vout, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vin, &iflag_);
    ASSERT_EQ(0,iflag_);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_same_vec_views_0_1__2_3)
  {
    if( !typeImplemented_ || problemTooSmall_ )
      return;

    // random data
    _ST_ shifts[2];
    for(int i = 0; i < 2; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0,iflag_);

    // safe vec1_
    SUBR(mvec_add_mvec)(st::one(), vec1_, st::zero(), vec2_, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(), vec1_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0,iflag_);
    // create views
    TYPE(mvec_ptr) vin = NULL;
    SUBR(mvec_view_block)(vec1_, &vin, 2, 3, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vout = NULL;
    SUBR(mvec_view_block)(vec1_, &vout, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vref = NULL;
    SUBR(mvec_view_block)(vec2_, &vref, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);

    // first generate reference data (safe calculation)
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A_, shifts, vin, beta, vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    // calculation (unsafe aliasing!)
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A_, shifts, vin, beta, vout, &iflag_);
    ASSERT_EQ(0,iflag_);

    // check vin
    SUBR(mvec_from_device)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_from_device)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_+VIDX(0,2,lda_),vec3_vp_+VIDX(0,2,lda_),nloc_,2,lda_,stride_,vflag_), mt::eps());
    // check result
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_,vec2_vp_,nloc_,2,lda_,stride_,vflag_), sqrt(mt::eps()));

    // delete views
    SUBR(mvec_delete)(vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vout, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vin, &iflag_);
    ASSERT_EQ(0,iflag_);
  }
#endif
#endif


#if MATNAME == MATNAME_sprandn
TEST_F(CLASSNAME,mvecT_times_mvec_after_spmvm)
{
  if( !typeImplemented_ || problemTooSmall_ )
    return;

  // fill vectors with ones
  SUBR(mvec_put_value)(vec1_,st::one(),&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_put_value)(vec2_,st::one(),&iflag_);
  ASSERT_EQ(0,iflag_);

  // perform an matrix-vector multiplication (to possibly fill up halo/padding with data)
  SUBR(sparseMat_times_mvec)(st::one(),A_,vec1_,st::one(),vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(sparseMat_times_mvec_communicate)(A_,vec2_,&iflag_);

  // fill columns with row/i and row/(i+1)
  // exploits sum_(k=1)^n 1/(k*(k+1)) = 1 - 1/(n+1)
  phist_gidx ilower;     
  phist_map_get_ilower(map_,&ilower,&iflag_);
  EXPECT_EQ(0,iflag_);// make sure this is a linear map, otherwise the
                      // analytical sum below is not matched
  for (int ii=0; ii< nloc_; ii++)
  {
    for (int j=0; j<nvec_; j++)
      vec2_vp_[VIDX(ii,j,lda_)] = (_ST_) -st::one()*(_MT_)((j+1)*1.0l/(ilower+ii+1));
    for (int i=0; i<nvec_; i++)
      vec1_vp_[VIDX(ii,i,lda_)] = (_ST_) st::one()*(_MT_)((i+1)*1.0l/(ilower+ii+2));
  }
  SUBR(mvec_to_device)(vec1_,&iflag_);
  SUBR(mvec_to_device)(vec2_,&iflag_);
  SUBR(mvecT_times_mvec)(st::one(),vec1_,vec2_,st::zero(),mat1_,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(sdMat_from_device)(mat1_,&iflag_);
  // check result
  sdMat_parallel_check(mat1_,&iflag_);
  ASSERT_EQ(0,iflag_);
  for(int i = 0; i < nvec_; i++)
  {
    for(int j = 0; j < nvec_; j++)
    {
      _MT_ val = -mt::one()*(i+1)*(j+1)*(1.0l - 1.0l/(_N_+1));
      ASSERT_NEAR(val, st::real(mat1_vp_[MIDX(i,j,m_lda_)]), _N_*20*mt::eps());
      ASSERT_NEAR(mt::zero(), st::imag(mat1_vp_[MIDX(i,j,m_lda_)]), _N_*20*mt::eps());
    }
  }
}
#endif


#if MATNAME==MATNAME_BENCH3D_8_A1

TEST_F(CLASSNAME,compare_with_rowFunc)
{
  if (MatNameEnumIsMatFunc(MATNAME)==false) return;
  if (!typeImplemented_ || problemTooSmall_) return;
  
  // arguments passed to the mvec* functions for filling the vectors
  int v_arg[2];
  v_arg[0]=_N_;
  v_arg[1]=_M_;

  // initialize vec1 and vec2 with row funcs
  SUBR(mvec_put_func)(vec1_,&PHIST_TG_PREFIX(mvec123func),v_arg,&iflag_);
  ASSERT_EQ(0,iflag_);
  SUBR(mvec_put_func)(vec2_,&PHIST_TG_PREFIX(mvec321func),v_arg,&iflag_);
  ASSERT_EQ(0,iflag_);

  _ST_ alpha = st::prand();
  _ST_ beta  = st::prand();
  SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  
  SUBR(mvec_from_device)(vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);
  
  // now construct the actual result in tmp_vec1 using only the defining functions.
  // we construct an additional vector here to make sure it uses the standard map of
  // the kernel library (typicaly a linear index space without load balancing or matrix-based partitioning)
  phist_lidx nloc = 0;
  phist_map_get_local_length(defaultMap_,&nloc,&iflag_);
      ASSERT_EQ(0,iflag_);

  TYPE(mvec_ptr) tmp_vec1=NULL;
  iflag_=0;
  // since we use temporary vector, we will only EXPECT instead of ASSERT so that
  // the memory is freed at the end of the function even if an assertion fails.
  PHISTTEST_MVEC_CREATE(&tmp_vec1, defaultMap_, nvec_, &iflag_);
  EXPECT_EQ(0,iflag_);
  _ST_* vec1_vp=NULL;
  phist_lidx lda_tmp;
  SUBR(mvec_extract_view)(tmp_vec1,&vec1_vp,&lda_tmp,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_gidx ilower;
  phist_map_get_ilower(defaultMap_,&ilower,&iflag_);
  EXPECT_EQ(0,iflag_);  
  bool row_func_error_encountered=false;
  for (int i=0; i<nloc_; i++)
  {
    phist_gidx row=ilower+i;
    for (int j=0; j<nvec_;j++)
    {
      _ST_ v2val;
      PHIST_TG_PREFIX(mvec321func)(row,j,&v2val,v_arg);
      vec1_vp[VIDX(i,j,lda_tmp)]=beta*v2val;
    }
    ghost_lidx len;
    ghost_gidx cols[_N_];
    ST val[_N_];
    iflag_=MATPDE3D_rowFunc(row,&len,cols,val,NULL);
    if (iflag_) {row_func_error_encountered=true; break;}
    for (int j=0; j<len; j++)
    {
      for (int k=0; k<nvec_; k++)
      {
        _ST_ v1val;
        iflag_=PHIST_TG_PREFIX(mvec123func)(cols[j],k,&v1val,v_arg);
        if (iflag_) {row_func_error_encountered=true; break;}
        vec1_vp[VIDX(i,k,lda_tmp)]+=alpha*val[j]*v1val;
      }
      if (row_func_error_encountered) break;
    }
    if (row_func_error_encountered) break;
  }

  SUBR(mvec_to_device)(tmp_vec1,&iflag_);
  EXPECT_EQ(0,iflag_);
  EXPECT_FALSE(row_func_error_encountered); // something wrong with row functions?

  // our expected result in tmp_vec1 is in the linear map defaultMap_, if necessary convert it to the permuted map map_
  SUBR(mvec_to_mvec)(tmp_vec1,vec1_,&iflag_);
  EXPECT_EQ(0,iflag_);
  
  EXPECT_NEAR(1.0,MvecsEqual(vec1_,vec2_),sqrt(mt::eps()));
  if (tmp_vec1!=NULL) 
  {
    SUBR(mvec_delete)(tmp_vec1,&iflag_);
    EXPECT_EQ(0,iflag_);
  }
}

TEST_F(CLASSNAME,fromRowFuncAndContext)
{
  if( !typeImplemented_ || problemTooSmall_ )
    return;

  TYPE(sparseMat_ptr) A1=NULL, A2=NULL;
  iflag_=PHIST_SPARSEMAT_QUIET;
  SUBR(sparseMat_create_fromRowFuncAndContext)(&A1,defaultContext_,7,&MATPDE3D_rowFunc,NULL,&iflag_);
  ASSERT_EQ(0,iflag_);
    iflag_=PHIST_SPARSEMAT_QUIET;
  SUBR(sparseMat_create_fromRowFuncAndContext)(&A2,context_,7,&MATPDE3D_rowFunc,NULL,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_const_map_ptr domain1, domain2, range1, range2, row1, row2, col1, col2;
  // first check if the created matrices are compatible with vectors of the given map:
  SUBR(sparseMat_get_range_map)(A1,&range1,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_domain_map)(A1,&domain1,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_maps_compatible(domain1,defaultMap_,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_maps_compatible(range1,defaultMap_,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_range_map)(A2,&range2,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_domain_map)(A2,&domain2,&iflag_);
  EXPECT_EQ(0,iflag_);
  
  phist_maps_compatible(domain2,map_,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_maps_compatible(range2,map_,&iflag_);
  EXPECT_EQ(0,iflag_);
  // then check if all of the maps are identical to those of A_ if the context of A is given
  SUBR(sparseMat_get_range_map)(A_,&range1,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_domain_map)(A_,&domain1,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_row_map)(A_,&row1,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_col_map)(A_,&col1,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_row_map)(A2,&row2,&iflag_);
  EXPECT_EQ(0,iflag_);
  SUBR(sparseMat_get_col_map)(A2,&col2,&iflag_);
  EXPECT_EQ(0,iflag_);
  //
  phist_maps_compatible(domain2,domain1,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_maps_compatible(range2,range1,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_maps_compatible(col2,col1,&iflag_);
  EXPECT_EQ(0,iflag_);
  phist_maps_compatible(row2,row1,&iflag_);
  EXPECT_EQ(0,iflag_);
  // delete matrices
  SUBR(sparseMat_delete)(A1,&iflag_);
  EXPECT_EQ(0,iflag_); 
  SUBR(sparseMat_delete)(A2,&iflag_);
  EXPECT_EQ(0,iflag_); 
}

#endif


#endif // DONT_INSTANTIATE

#ifdef DONT_INSTANTIATE
#undef DONT_INSTANTIATE
#endif
