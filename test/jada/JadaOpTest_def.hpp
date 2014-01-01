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

      if (typeImplemented_)
      {
        SUBR(mvec_create)(&q_,map_,_NVP_,&ierr_);
        ASSERT_EQ(0,ierr_);
        sigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
          sigma_[i] = st::rand();

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&ierr_);
        ASSERT_EQ(0,ierr_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_QR)(q_,Rtmp,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(sdMat_delete)(Rtmp,&ierr_);
        ASSERT_EQ(0,ierr_);

        SUBR(read_mat)("sprandn",nglob_,&A1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        opA1_ = new TYPE(op);
        SUBR(op_wrap_crsMat)(opA1_, A1_, &ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(read_mat)("speye",nglob_,&I_,&ierr_);
        ASSERT_EQ(0,ierr_);
        opI_ = new TYPE(op);
        SUBR(op_wrap_crsMat)(opI_, I_, &ierr_);
        ASSERT_EQ(0,ierr_);
      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      MTest::TearDown();
      VTest::TearDown();
      if (typeImplemented_)
      {
        delete opI_;
        SUBR(crsMat_delete)(A1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        delete opA1_;
        SUBR(crsMat_delete)(I_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(q_,&ierr_);
        ASSERT_EQ(0,ierr_);
        delete[] sigma_;
      }
    }

    TYPE(crsMat_ptr) A1_; 
    TYPE(op_ptr) opA1_;
    TYPE(crsMat_ptr) I_; 
    TYPE(op_ptr) opI_;
    TYPE(mvec_ptr) q_;
    _ST_* sigma_;
};


  TEST_F(CLASSNAME, create_and_delete)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, create_and_delete_generalized)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,opI_,q_,q_,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, apply_only_projection)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opI_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_random)(vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      // apply
      for(int i = 0; i < _NV_; i++)
        sigma_[i] = st::zero();
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

      {
        TYPE(mvec_ptr) AX = NULL;
        SUBR(jadaOp_view_AX)(jdOp.A,&AX,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_add_mvec)(st::one(),AX,st::zero(),vec1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(AX,&ierr_);
        ASSERT_EQ(0,ierr_);
      }

      // vec1_ = I*vec2_ ?
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif

      // vec3_ = (I-q_*q_') vec2_ ?
      SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(-st::one(),vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()));
#endif

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, apply_shifted_projection)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opI_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_random)(vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

      {
        TYPE(mvec_ptr) AX = NULL;
        SUBR(jadaOp_view_AX)(jdOp.A,&AX,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_add_mvec)(st::one(),AX,st::zero(),vec1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(AX,&ierr_);
        ASSERT_EQ(0,ierr_);
      }

      // vec1_ = I*vec2_ ?
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif

      // vec3_ = (I-q_*q_') (I+(sigma_i))vec2_ ?
      for(int i = 0; i < _NV_; i++)
        sigma_[i] += st::one();
      SUBR(mvec_vscale)(vec2_,sigma_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvecT_times_mvec)(st::one(),q_,vec2_,st::zero(),mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(-st::one(),vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()));
#endif

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, apply_check_AX)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_random)(vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

      {
        TYPE(mvec_ptr) AX = NULL;
        SUBR(jadaOp_view_AX)(jdOp.A,&AX,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_add_mvec)(st::one(),AX,st::zero(),vec1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(AX,&ierr_);
        ASSERT_EQ(0,ierr_);
      }

      // vec1_ = A1_*vec2_ ?
      SUBR(crsMat_times_mvec)(st::one(),A1_,vec2_,st::zero(),vec3_,&ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec3_vp_,vec1_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec3_vp_,vec1_vp_,nloc_,nvec_,lda_,stride_));
#endif

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, apply_check_result)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_random)(vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec2_,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

      {
        TYPE(mvec_ptr) AX = NULL;
        SUBR(jadaOp_view_AX)(jdOp.A,&AX,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_add_mvec)(st::one(),AX,st::zero(),vec1_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(AX,&ierr_);
        ASSERT_EQ(0,ierr_);
      }

      // vec3_ = (I-q_*q_') (vec1+sigma_i*I)vec2_ ?
      SUBR(mvec_vadd_mvec)(sigma_,vec2_,st::one(),vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvecT_times_mvec)(st::one(),q_,vec1_,st::zero(),mat2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(-st::one(),vec1_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_times_sdMat)(st::one(),q_,mat2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, apply_test_alpha)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      TYPE(mvec_ptr) vec4;
      SUBR(mvec_create)(&vec4,map_,_NV_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_random)(vec4,&ierr_);
      ASSERT_EQ(0,ierr_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec4,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

      TYPE(mvec_ptr) vec5;
      SUBR(mvec_create)(&vec5,map_,_NV_,&ierr_);
      ASSERT_EQ(0,ierr_);
      _ST_ alpha = st::rand();
      jdOp.apply(alpha,jdOp.A,vec4,st::zero(),vec5,&ierr_);
      ASSERT_EQ(0,ierr_);

      // vec5 = alpha*vec3_
      SUBR(mvec_add_mvec)(-st::one(),vec5,alpha,vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(mvec_delete)(vec5,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(vec4,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }


  TEST_F(CLASSNAME, apply_test_beta)
  {
    if( typeImplemented_ )
    {
      TYPE(op) jdOp;
      SUBR(jadaOp_create)(opA1_,NULL,q_,NULL,sigma_,_NV_,&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);

      TYPE(mvec_ptr) vec4;
      SUBR(mvec_create)(&vec4,map_,_NV_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_random)(vec4,&ierr_);
      ASSERT_EQ(0,ierr_);

      // apply
      jdOp.apply(st::one(),jdOp.A,vec4,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

      TYPE(mvec_ptr) vec5;
      SUBR(mvec_create)(&vec5,map_,_NV_,&ierr_);
      SUBR(mvec_put_value)(vec5,st::one(),&ierr_);
      ASSERT_EQ(0,ierr_);
      // we assume vec5 in q^orth, so make it so:
      SUBR(mvecT_times_mvec)(st::one(),q_,vec5,st::zero(),mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_times_sdMat)(-st::one(),q_,mat1_,st::one(),vec5,&ierr_);
      ASSERT_EQ(0,ierr_);
      // add beta (I-qq')*ONE to vec3_
      _ST_ beta = st::rand();
      SUBR(mvec_add_mvec)(beta,vec5,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      jdOp.apply(st::one(),jdOp.A,vec4,beta,vec5,&ierr_);
      ASSERT_EQ(0,ierr_);

      // vec5 == vec3_?
      SUBR(mvec_add_mvec)(st::one(),vec5,-st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);

#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nvec_,nloc_,lda_,stride_,st::zero()),10*VTest::releps());
#else
      ASSERT_NEAR(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero()),10*VTest::releps());
#endif

      SUBR(mvec_delete)(vec5,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(vec4,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(jadaOp_delete)(&jdOp,&ierr_);
      ASSERT_EQ(0,ierr_);
    }
  }

