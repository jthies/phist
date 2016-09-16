#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithMassMat<_ST_,_N_>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,3>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>
{

  public:
    typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
    typedef KernelTestWithMassMat<_ST_,_N_> BTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,3> VTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;

    static void SetUpTestCase()
    {
      int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
      SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
      BTest::SetUpTestCase(map_);
      VTest::SetUpTestCase();
      MTest::SetUpTestCase();
    }

    /*! Set up routine.
    */
    virtual void SetUp()
    {
      SparseMatTest::SetUp();
      BTest::SetUp();
      VTest::SetUp();
      MTest::SetUp();

      if( typeImplemented_ && !problemTooSmall_ )
      {

        // disable the test because TSQR will not work.
        // This situation occurs if we have a small matrix (_N_=25, say)
        // and many Q vectors (e.g. 10) with multiple MPI procs.
        int globalTooSmall = _N_ < _NVP_;
#ifdef PHIST_HAVE_MPI
        int localTooSmall = nloc_ < _NVP_;
        iflag_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        ASSERT_EQ(0,iflag_);
#endif
        problemTooSmall_ = globalTooSmall != 0;
      }

      if (typeImplemented_ && !problemTooSmall_)
      {

        SUBR(mvec_random)(vec1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_random)(vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);

        PHISTTEST_MVEC_CREATE(&q_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        PHISTTEST_MVEC_CREATE(&qb_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        PHISTTEST_MVEC_CREATE(&Bq_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        sigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
          sigma_[i] = st::prand();

        opA_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opA_, A_, &iflag_);
        ASSERT_EQ(0,iflag_);

        opB_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opB_, B_, &iflag_);
        ASSERT_EQ(0,iflag_);

        opAB_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat_pair)(opAB_, A_, B_, &iflag_);
        ASSERT_EQ(0,iflag_);

        // create random orthogonal Q and B-orhtogonalQ_B
        SUBR(mvec_random)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_random)(qb_,&iflag_);
        ASSERT_EQ(0,iflag_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        int rankQ=0;
        SUBR(orthog)(NULL,q_,NULL,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(orthog)(NULL,qb_,opB_,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(sdMat_delete)(Rtmp,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        SUBR(sparseMat_times_mvec)(st::one(),B_,qb_,st::zero(),Bq_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      if (typeImplemented_ && !problemTooSmall_)
      {
        if( opA_ != NULL )
        {
          opA_->destroy(opA_,&iflag_);
          ASSERT_EQ(0,iflag_);
          delete opA_;
        }
        opA_ = NULL;
        if( opB_ != NULL )
        {
          opB_->destroy(opB_,&iflag_);
          ASSERT_EQ(0,iflag_);
          delete opB_;
        }
        opB_ = NULL;
        if( opAB_ != NULL )
        {
          opAB_->destroy(opAB_,&iflag_);
          ASSERT_EQ(0,iflag_);
          delete opAB_;
        }
        opA_ = NULL;
        SUBR(mvec_delete)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(qb_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(Bq_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( sigma_ != NULL)
          delete[] sigma_;
        sigma_ = NULL;
      }

      MTest::TearDown();
      VTest::TearDown();
      SparseMatTest::TearDown();
      BTest::TearDown();
    }

    static void TearDownTestCase()
    {
      MTest::TearDownTestCase();
      VTest::TearDownTestCase();
      SparseMatTest::TearDownTestCase();
      BTest::TearDownTestCase();
    }

    TYPE(linearOp_ptr) opA_ = NULL, opB_ = NULL, opAB_ = NULL;
    TYPE(mvec_ptr) q_ = NULL, Bq_ = NULL, qb_ = NULL;
    _ST_* sigma_ = NULL;
};


  TEST_F(CLASSNAME, create_and_delete)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  TEST_F(CLASSNAME, create_and_delete_withB)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opAB_,opB_,qb_,Bq_,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


#if MATNAME == MATNAME_speye
  TEST_F(CLASSNAME, apply_only_projection)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      for(int i = 0; i < _NV_; i++) sigma_[i] = st::zero();
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);


      // vec3_ = (I-q_*q_') vec2_ ?
      SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec3_,st::zero()));

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  TEST_F(CLASSNAME, apply_only_projection_withB)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      for(int i = 0; i < _NV_; i++) sigma_[i] = st::zero();
      SUBR(jadaOp_create)(opAB_,opB_,qb_,Bq_,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
      // apply, vec1 = (I-Bqq') I (I-qBq')vec2
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // vec3_ = (I-Bq_*q_') vec2_. Note that vec3=vec2 by construction
      SUBR(mvecT_times_mvec)(st::one(),Bq_,vec3_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(-st::one(),qb_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // vec1_ <- vec1 - (vec2 projected twice) should be zero
      SUBR(mvecT_times_mvec)(st::one(),qb_,vec3_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec3_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),Bq_,mat2_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec1_,st::zero()));

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_shifted_projection)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);


      // vec3_ = (I-q_*q_') (I+(sigma_i))vec2_ ?
      for(int i = 0; i < _NV_; i++)
        sigma_[i] += st::one();
      SUBR(mvec_vscale)(vec2_,sigma_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),MvecEqual(vec3_,st::zero()),1000*mt::eps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
#endif


  TEST_F(CLASSNAME, apply_check_result)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);


      // vec3_ = (I-q_*q_') (vec1+sigma_i*I)vec2_ ?
      opA_->apply(st::one(), opA_->A, vec2_, st::zero(), vec1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_vadd_mvec)(sigma_,vec2_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),q_,vec1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec1_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MvecEqual(vec3_,st::zero()),10*VTest::releps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  TEST_F(CLASSNAME, apply_check_result_withB)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opAB_,opB_,qb_,Bq_,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply, vec2 = (jdOp)*vec3, save vec2 for the final comparison
      jdOp.apply(st::one(),jdOp.A,vec3_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);


      // vec3_ = (I-q_*Bq_') vec3_.
      SUBR(mvecT_times_mvec)(st::one(),Bq_,vec3_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(-st::one(),qb_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // vec1_ = (A+sigma_i*B)vec3
      SUBR(sparseMat_times_mvec)(st::one(),B_,vec3_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_vscale)(vec1_,sigma_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sparseMat_times_mvec)(st::one(), A_, vec3_, +st::one(), vec1_, &iflag_);
      ASSERT_EQ(0,iflag_);

      // vec1 = (I - Bqq')vec1
      SUBR(mvecT_times_mvec)(st::one(),qb_,vec1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(-st::one(),Bq_,mat2_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MvecsEqual(vec1_,vec2_),100.0*VTest::releps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_test_alpha)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec4;
      PHISTTEST_MVEC_CREATE(&vec4,map_,_NV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec4,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec5;
      PHISTTEST_MVEC_CREATE(&vec5,map_,_NV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = st::prand();
      jdOp.apply(alpha,jdOp.A,vec4,st::zero(),vec5,&iflag_);
      ASSERT_EQ(0,iflag_);

      // vec5 = alpha*vec3_
      SUBR(mvec_add_mvec)(-st::one(),vec5,alpha,vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MvecEqual(vec3_,st::zero()),10*VTest::releps());

      SUBR(mvec_delete)(vec5,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_test_beta)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec4;
      PHISTTEST_MVEC_CREATE(&vec4,map_,_NV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec4,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec5;
      PHISTTEST_MVEC_CREATE(&vec5,map_,_NV_,&iflag_);
      SUBR(mvec_put_value)(vec5,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      // we assume vec5 in q^orth, so make it so:
      SUBR(mvecT_times_mvec)(st::one(),q_,vec5,st::zero(),mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(-st::one(),q_,mat1_,st::one(),vec5,&iflag_);
      ASSERT_EQ(0,iflag_);
      // add beta (I-qq')*ONE to vec3_
      _ST_ beta = st::prand();
      SUBR(mvec_add_mvec)(beta,vec5,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      jdOp.apply(st::one(),jdOp.A,vec4,beta,vec5,&iflag_);
      ASSERT_EQ(0,iflag_);

      // vec5 == vec3_?
      SUBR(mvec_add_mvec)(st::one(),vec5,-st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MvecEqual(vec3_,st::zero()),10*VTest::releps());

      SUBR(mvec_delete)(vec5,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  // checks that A_orth Q = 0 
  // NOTE: test disabled because this only holds if we pre-project,
  // right now our jadaOp is defined as (I-VV')(A-sigma*I) if B=I, however.
  TEST_F(CLASSNAME, DISABLED_orthogonality_of_jadaOp)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      _ST_ sigma[_NVP_];
      for (int i=0; i<_NVP_; i++) sigma[i]=st::prand();
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma,_NVP_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply to Q itself (should give 0)
      jdOp.apply(st::one(),jdOp.A,q_,st::zero(),qb_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(st::one(),MvecEqual(qb_,st::zero()),VTest::releps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  // checks that A_orth Q = 0 and Q'A_orth X = 0
  TEST_F(CLASSNAME, check_orthogonality_of_output_vector)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opA_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply to X
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check that result is orthogonal to q
      SUBR(mvecT_times_mvec)(st::one(),q_,vec3_,st::zero(),mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(st::one(),SdMatEqual(mat1_,st::zero()),15*mt::eps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  // checks that A_orth*Q = 0
  TEST_F(CLASSNAME, check_orthogonality_of_jadaOp_withB)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      _ST_ sigma[_NVP_];
      for (int i=0; i<_NVP_; i++) sigma[i]=st::prand();
      SUBR(jadaOp_create)(opAB_,opB_,qb_,Bq_,sigma,_NVP_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply to BQ itself (should give 0)
      jdOp.apply(st::one(),jdOp.A,qb_,st::zero(),q_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(st::one(),MvecEqual(q_,st::zero()),10*VTest::releps());

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check that result is orthogonal to q
      SUBR(mvecT_times_mvec)(st::one(),qb_,vec3_,st::zero(),mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(st::one(),SdMatEqual(mat1_,st::zero()),20*mt::eps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  TEST_F(CLASSNAME, check_orthogonality_of_output_vector_withB)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(linearOp) jdOp;
      SUBR(jadaOp_create)(opAB_,opB_,qb_,Bq_,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply to X
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // check that result is orthogonal to q
      SUBR(mvecT_times_mvec)(st::one(),qb_,vec3_,st::zero(),mat1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_NEAR(st::one(),SdMatEqual(mat1_,st::zero()),10*mt::eps());

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }
