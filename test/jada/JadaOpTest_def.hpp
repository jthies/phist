#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>
{

  public:
    typedef KernelTestWithVectors<_ST_,_N_,_NV_> VTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NVP_> VProjTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;


    /*! Set up routine.
    */
    virtual void SetUp()
    {
      VTest::SetUp();
      MTest::SetUp();

      if( typeImplemented_ )
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
        typeImplemented_ = typeImplemented_ && !globalTooSmall;
      }

      if (typeImplemented_)
      {
        SUBR(mvec_create)(&q_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        sigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
          sigma_[i] = st::prand();

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        int rankQ=0;
        SUBR(orthog)(NULL,q_,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(sdMat_delete)(Rtmp,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(read_mat)("sprandn",nglob_,&A1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        opA1_ = new TYPE(op);
        SUBR(op_wrap_crsMat)(opA1_, A1_, &iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(read_mat)("speye",nglob_,&I_,&iflag_);
        ASSERT_EQ(0,iflag_);
        opI_ = new TYPE(op);
        SUBR(op_wrap_crsMat)(opI_, I_, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    virtual void replaceMap(const_map_ptr_t map)
    {
      if (typeImplemented_)
      {
        // delete old vecs
        SUBR(mvec_delete)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);

        VTest::replaceMap(map);

        // create new vecs

        SUBR(mvec_create)(&q_,map_,_NVP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        sigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
          sigma_[i] = st::prand();

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        int rankQ=0;
        SUBR(orthog)(NULL,q_,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(sdMat_delete)(Rtmp,&iflag_);
        ASSERT_EQ(0,iflag_);

      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      if (typeImplemented_)
      {
        if( opI_ != NULL )
          delete opI_;
        SUBR(crsMat_delete)(A1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( opA1_ != NULL )
          delete opA1_;
        SUBR(crsMat_delete)(I_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( sigma_ != NULL)
          delete[] sigma_;
      }

      MTest::TearDown();
      VTest::TearDown();
    }

    TYPE(crsMat_ptr) A1_ = NULL;
    TYPE(op_ptr) opA1_ = NULL;
    TYPE(crsMat_ptr) I_ = NULL;
    TYPE(op_ptr) opI_ = NULL;
    TYPE(mvec_ptr) q_ = NULL;
    _ST_* sigma_ = NULL;
};


  TEST_F(CLASSNAME, create_and_delete)
  {
    if( typeImplemented_ )
    {
      // replace map to by the map of the current matrix
      const_map_ptr_t map = NULL;
      SUBR(crsMat_get_domain_map)(A1_,&map,&iflag_);
      ASSERT_EQ(0,iflag_);
      replaceMap(map);
      ASSERT_EQ(0,iflag_);

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, DISABLE_create_and_delete_generalized)
  {
    if( typeImplemented_ )
    {
      // we need fitting maps??

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,opI_,q_,q_,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_only_projection)
  {
    if( typeImplemented_ )
    {
      // replace map to by the map of the current matrix
      const_map_ptr_t map = NULL;
      SUBR(crsMat_get_domain_map)(I_,&map,&iflag_);
      ASSERT_EQ(0,iflag_);
      replaceMap(map);
      ASSERT_EQ(0,iflag_);

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opI_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      for(int i = 0; i < _NV_; i++)
        sigma_[i] = st::zero();
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);


      // vec3_ = (I-q_*q_') vec2_ ?
      SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()));
#endif

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_shifted_projection)
  {
    if( typeImplemented_ )
    {
      // replace map to by the map of the current matrix
      const_map_ptr_t map = NULL;
      SUBR(crsMat_get_domain_map)(I_,&map,&iflag_);
      ASSERT_EQ(0,iflag_);
      replaceMap(map);
      ASSERT_EQ(0,iflag_);

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opI_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_random)(vec2_,&iflag_);
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
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()));
#endif

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_check_result)
  {
    if( typeImplemented_ )
    {
      // replace map to by the map of the current matrix
      const_map_ptr_t map = NULL;
      SUBR(crsMat_get_domain_map)(A1_,&map,&iflag_);
      ASSERT_EQ(0,iflag_);
      replaceMap(map);
      ASSERT_EQ(0,iflag_);

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);


      // vec3_ = (I-q_*q_') (vec1+sigma_i*I)vec2_ ?
      opA1_->apply(st::one(), opA1_->A, vec2_, st::zero(), vec1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_vadd_mvec)(sigma_,vec2_,st::one(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),q_,vec1_,st::zero(),mat2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec1_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }


  TEST_F(CLASSNAME, apply_test_alpha)
  {
    if( typeImplemented_ )
    {
      // replace map to by the map of the current matrix
      const_map_ptr_t map = NULL;
      SUBR(crsMat_get_domain_map)(A1_,&map,&iflag_);
      ASSERT_EQ(0,iflag_);
      replaceMap(map);
      ASSERT_EQ(0,iflag_);

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec4;
      SUBR(mvec_create)(&vec4,map_,_NV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec4,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec5;
      SUBR(mvec_create)(&vec5,map_,_NV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ alpha = st::prand();
      jdOp.apply(alpha,jdOp.A,vec4,st::zero(),vec5,&iflag_);
      ASSERT_EQ(0,iflag_);

      // vec5 = alpha*vec3_
      SUBR(mvec_add_mvec)(-st::one(),vec5,alpha,vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

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
    if( typeImplemented_ )
    {
      // replace map to by the map of the current matrix
      const_map_ptr_t map = NULL;
      SUBR(crsMat_get_domain_map)(A1_,&map,&iflag_);
      ASSERT_EQ(0,iflag_);
      replaceMap(map);
      ASSERT_EQ(0,iflag_);

      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec4;
      SUBR(mvec_create)(&vec4,map_,_NV_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec4,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      TYPE(mvec_ptr) vec5;
      SUBR(mvec_create)(&vec5,map_,_NV_,&iflag_);
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

#ifdef PHIST_MVECS_ROW_MAJOR
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(mvec_delete)(vec5,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_delete)(vec4,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(jadaOp_delete)(&jdOp,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

