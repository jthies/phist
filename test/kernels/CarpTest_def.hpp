#include "../tools/TestHelpers.h"

#ifndef CLASSNAME
#error "file not included correctly."
#endif

#if !defined(IS_COMPLEX)&&defined(IS_DOUBLE)

/*! Test fixure. 
  
  basic tests for CARP kernel, creates a test matrix defined by
  _MATNAME_ and some vectors, checks if the kernel
  works with and without RHS, with and without real or complex shift.
  
  It is not easy to rigorously test this kernel routine because we give
  the kernel lib quite some freedom in deciding e.g. on the order of Kaczmarz
  updates and the parallelization scheme. What we can test is
  
  * every vector element should be updated
  * operator should be symmetric (hermitian for complex shift), S=X'(DKSWP(X))=S'
  
  given an Identity matrix,
  
  * for rhs 0 and shift 0  x=>0
  * for rhs 0 and shift 1  x=>x
  ... TODO: add tests ...
  
*/
class CLASSNAME: public KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public KernelTestWithVectors<_ST_,_N_,_NV_,0,3> 
{

public:
  typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,3>  VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_> MTest;

  typedef TestWithType<MT> MT_Test;

  static void SetUpTestCase()
  {
    int sparseMatCreateFlag=PHIST_SPARSEMAT_OPT_CARP | getSparseMatCreateFlag(_N_,_NV_);
    SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
    VTest::SetUpTestCase();
    MT_Test::SetUpTestCase();
    
  }

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    SparseMatTest::SetUp();
    VTest::SetUp();

    
    vec1b_=NULL;
    PHISTTEST_MVEC_CREATE(&vec1b_,map_,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    vec2b_=NULL; 
    PHISTTEST_MVEC_CREATE(&vec2b_,map_,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    vec3b_=NULL;
    PHISTTEST_MVEC_CREATE(&vec3b_,map_,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    // set pointers for the tests
    x_r=vec1_; x_r_bak=vec1b_;
    x_i=vec2_; x_i_bak=vec2b_;
    b=vec3_;

    I_=NULL;
    if (typeImplemented_ && !problemTooSmall_)
    {
      iflag_=PHIST_SPARSEMAT_OPT_CARP | getSparseMatCreateFlag(_N_,_NV_);
      SUBR(sparseMat_create_fromRowFunc)(&I_,comm_,_N_,_N_,1,&SUBR(idfunc),NULL,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec1b_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec2b_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec3b_,&iflag_);
      ASSERT_EQ(0,iflag_);

      iflag_=PHIST_MVEC_REPLICATE_DEVICE_MEM;
      x_vec1_=new TYPE(x_mvec)(vec1_,vec1b_,0,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=PHIST_MVEC_REPLICATE_DEVICE_MEM;
      x_vec2_=new TYPE(x_mvec)(vec2_,vec2b_,0,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=PHIST_MVEC_REPLICATE_DEVICE_MEM;
      x_vec3_=new TYPE(x_mvec)(vec3_,vec3b_,0,&iflag_);
      ASSERT_EQ(0,iflag_);


      iflag_=PHIST_MVEC_REPLICATE_DEVICE_MEM;
      phist_Zmvec_create(&z_vec1_,map_,nvec_,&iflag_);
      cTypeImplemented_=(iflag_!=PHIST_NOT_IMPLEMENTED);         

      sigma_[0]=1.0-ct::cmplx_I();        minus_sigma_[0]=-sigma_[0];
      sigma_r_[0]=ct::real(sigma_[0]);     sigma_i_[0]=ct::imag(sigma_[0]);
      for (int i=1; i<nvec_; i++)
      {
        sigma_[i]=ct::prand();                 minus_sigma_[i]=-sigma_[i];
        sigma_r_[i]=ct::real(sigma_[i]);      sigma_i_[i]=ct::imag(sigma_[i]);
      }

      for (int i=0; i<nvec_; i++) omega_[i]=1.84299;



      if (cTypeImplemented_)
      {
        ASSERT_EQ(0,iflag_);
        iflag_=PHIST_MVEC_REPLICATE_DEVICE_MEM;
        phist_Zmvec_create(&z_vec2_,map_,nvec_,&iflag_); \
        ASSERT_EQ(0,iflag_);
        iflag_=PHIST_MVEC_REPLICATE_DEVICE_MEM;
        phist_Zmvec_create(&z_vec3_,map_,nvec_,&iflag_); \
        ASSERT_EQ(0,iflag_);
        MvecCopyX2Z(x_vec1_,z_vec1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        MvecCopyX2Z(x_vec2_,z_vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        MvecCopyX2Z(x_vec3_,z_vec3_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // note: we make sure the complex matrices use the same map as the real ones, this sames some
        // trouble when comparing result vectors.
        iflag_=PHIST_SPARSEMAT_QUIET;
        phist_ZsparseMat_create_fromRowFuncAndMap(&z_A_,map_,7,&ZMATFUNC,NULL,&iflag_);
        ASSERT_EQ(0,iflag_);
        iflag_=PHIST_SPARSEMAT_QUIET;
        phist_ZsparseMat_create_fromRowFuncAndMap(&z_A_shift0_,map_,7,&ZMATFUNC,&sigma_[0],&iflag_);
        ASSERT_EQ(0,iflag_);
        iflag_=PHIST_SPARSEMAT_QUIET;
        phist_ZsparseMat_create_fromRowFuncAndMap(&z_A_shift1_,map_,7,&ZMATFUNC,&sigma_[1],&iflag_);
        ASSERT_EQ(0,iflag_);

      }
      else
      {
        z_vec1_=NULL;
        z_vec2_=NULL;
        z_vec3_=NULL;
        z_A_=NULL;
        z_A_shift0_=NULL;
        z_A_shift1_=NULL;
      }
    
      // check if CARP is implemented at all:
      void* work;
      SUBR(carp_setup)(I_, 1, sigma_r_,
          &work, &iflag_);
      if (iflag_==PHIST_NOT_IMPLEMENTED) 
      {
        iflag_=0;
        carpImplemented_=false;
      }
      else 
      {
        carpImplemented_=true;
        iflag_=0;
        SUBR(carp_destroy)(I_, work, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
#ifndef IS_COMPLEX
      SUBR(carp_setup_rc)(I_, 1, sigma_r_,sigma_i_,
          &work, &iflag_);
      if (iflag_==PHIST_NOT_IMPLEMENTED) 
      {
        iflag_=0;
        rc_carpImplemented_=false;
      }
      else 
      {
        rc_carpImplemented_=true;
        iflag_=0;
        SUBR(carp_destroy)(I_, work, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
#endif
      x_A_=new TYPE(x_sparseMat);
      x_A_->A_=A_;
      x_A_->sigma_r_=new _ST_[nvec_];
      x_A_->sigma_i_=new _ST_[nvec_];
      x_A_->Vproj_=NULL;
      for (int i=0; i<nvec_;i++)
      {
        x_A_->sigma_r_[i]=sigma_r_[i];
        x_A_->sigma_i_[i]=sigma_i_[i];
      }
    }
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      if (vec1b_!=NULL)
      {
        SUBR(mvec_delete)(vec1b_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if (vec2b_!=NULL)
      {
        SUBR(mvec_delete)(vec2b_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if (vec3b_!=NULL)
      {
        SUBR(mvec_delete)(vec3b_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      ASSERT_EQ(0,delete_mat(I_));

      if (x_A_ != NULL)
      {
        delete [] x_A_->sigma_r_;
        delete [] x_A_->sigma_i_;
        delete x_A_;
      }
      
      if (x_vec1_!=NULL) delete x_vec1_;
      if (x_vec2_!=NULL) delete x_vec2_;
      if (x_vec3_!=NULL) delete x_vec3_;
  
      if (cTypeImplemented_)
      {
        if (z_vec1_!=NULL)
        {
          phist_Zmvec_delete(z_vec1_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        if (z_vec2_!=NULL)
        {
          phist_Zmvec_delete(z_vec2_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        if (z_vec3_!=NULL)
        {
          phist_Zmvec_delete(z_vec3_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        if (z_A_!=NULL)
        {
          phist_ZsparseMat_delete(z_A_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        if (z_A_shift0_!=NULL)
        {
          phist_ZsparseMat_delete(z_A_shift0_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        if (z_A_shift1_!=NULL)
        {
          phist_ZsparseMat_delete(z_A_shift1_,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
      }
    }
    
    VTest::TearDown();
    SparseMatTest::TearDown();
  }

  static void TearDownTestCase()
  {
    VTest::TearDownTestCase();
    SparseMatTest::TearDownTestCase();
  }



  // helper function to apply [x_r, x_i] = dkswp(A-sigma[j]*I, b, x_r, x_i)
  // and copy original X vectors to x_r_bak, x_i_bak.
  void create_and_apply_carp(TYPE(const_sparseMat_ptr) A)
  {

    SUBR(mvec_add_mvec)(st::one(),x_r,st::zero(),x_r_bak,&iflag_);
    ASSERT_EQ(0,iflag_);

    void* work;
    SUBR(carp_setup)(A, _NV_, sigma_r_, &work, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    SUBR(carp_sweep)(A, sigma_r_,b,x_r, work, omega_, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    SUBR(carp_destroy)(A, work, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    return;
  
  }

#ifndef IS_COMPLEX
  // helper function to apply [x_r, x_i] = dkswp(A-sigma[j]*I, b, x_r, x_i)
  // and copy original X vectors to x_r_bak, x_i_bak.
  void create_and_apply_carp_rc(TYPE(const_sparseMat_ptr) A)
  {

    SUBR(mvec_add_mvec)(st::one(),x_r,st::zero(),x_r_bak,&iflag_);
    ASSERT_EQ(0,iflag_);

    if (x_i!=NULL)
    {
      SUBR(mvec_add_mvec)(st::one(),x_i,st::zero(),x_i_bak,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    void* work;
    SUBR(carp_setup_rc)(A, _NV_, sigma_r_,sigma_i_, &work, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    SUBR(carp_sweep_rc)(A, sigma_r_,sigma_i_,b,x_r,x_i, work, omega_, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    SUBR(carp_destroy)(A, work, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    return;
  }
#endif
protected:

  int delete_mat(TYPE(sparseMat_ptr) &A)
  {
    if (A!=NULL)
    {
      SUBR(sparseMat_delete)(A,&iflag_);
      A = NULL;
    }
    return iflag_;
  }

  void check_symmetry(TYPE(const_mvec_ptr) X, TYPE(const_mvec_ptr) OPX,_MT_ tol=10*mt::eps())
  {
    _MT_ max_err=mt::zero();
    if (typeImplemented_ && !problemTooSmall_)
    {
      TYPE(sdMat_ptr) M=NULL;
      SUBR(sdMat_create)(&M,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),X,OPX,st::zero(),M,&iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ normF;
      SUBR(sdMat_normF)(M,&normF,&iflag_);
      EXPECT_FALSE(normF==mt::zero());
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(sdMat_from_device)(M,&iflag_);
      PHIST_SOUT(PHIST_DEBUG,"CARP: X'*OP*X=\n");
      SUBR(sdMat_print)(M,&iflag_);
#endif
      max_err=MTest::symmetry_check(M,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_delete)(M,&iflag_);
      ASSERT_EQ(0,iflag_);
      EXPECT_NEAR(mt::one(),max_err+mt::one(),tol);
    }
  }

  void check_symmetry_rc(TYPE(const_mvec_ptr) X_r, TYPE(const_mvec_ptr) X_i, 
                            TYPE(const_mvec_ptr) OPX_r, TYPE(const_mvec_ptr) OPX_i, _MT_ tol=10*mt::eps())
  {
    _MT_ max_err_r=mt::zero();
    _MT_ max_err_i=mt::zero();
    if (typeImplemented_ && !problemTooSmall_)
    {
      TYPE(sdMat_ptr) M_r=NULL, M_i=NULL;
      SUBR(sdMat_create)(&M_r,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_create)(&M_i,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(mvecT_times_mvec)(st::one(),X_r,OPX_r,st::zero(),M_r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(+st::one(),X_i,OPX_i,st::one(),M_r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),X_r,OPX_i,st::zero(),M_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(-st::one(),X_i,OPX_r,st::one(),M_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M_r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M_i,&iflag_);
      ASSERT_EQ(0,iflag_);
      phist_lidx ldmr,ldmi;
      _ST_ *M_r_raw=NULL, *M_i_raw=NULL;
      SUBR(sdMat_extract_view)(M_r,&M_r_raw,&ldmr,&iflag_);
      SUBR(sdMat_extract_view)(M_i,&M_i_raw,&ldmi,&iflag_);
      for (int i=0; i<_NV_; i++)
      {
        for (int j=i+1; j<_NV_;j++)
        {
          PHIST_SOUT(PHIST_INFO,"%d %d %25.16e%+25.16ei\n",i,j,M_r_raw[i*ldmr+j],M_i_raw[i*ldmi+j]);
          PHIST_SOUT(PHIST_INFO,"%d %d %25.16e%+25.16ei\n",j,i,M_r_raw[j*ldmr+i],M_i_raw[j*ldmi+i]);
          max_err_r = std::max(max_err_r,std::abs(M_r_raw[i*ldmr+j]-M_r_raw[j*ldmr+i]));
          max_err_i = std::max(max_err_i,std::abs(M_i_raw[i*ldmi+j]+M_i_raw[j*ldmi+i]));
        }
      }
      SUBR(sdMat_delete)(M_r,&iflag_);
      SUBR(sdMat_delete)(M_i,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    ASSERT_NEAR(mt::one(),max_err_r+mt::one(),tol);
    ASSERT_NEAR(mt::one(),max_err_i+mt::one(),tol);
  }

  void do_spmv_test(double alpha, double beta)
  {
    // sanity check of initial status
    ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
    ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_));
    
    phist_d_complex z_alpha = (phist_d_complex)alpha;
    phist_d_complex z_beta = (phist_d_complex)beta;
  
    for (int i=0; i<nvec_; i++)
    {
      minus_sigma_[i]=-(x_A_->sigma_r_[i] + ct::cmplx_I()*x_A_->sigma_i_[i]);
    }
  
    phist_ZsparseMat_times_mvec_vadd_mvec(z_alpha, z_A_, minus_sigma_, z_vec1_, z_beta, z_vec2_, &iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(x_sparseMat_times_mvec)(alpha, x_A_, x_vec1_, beta, x_vec2_, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
    ASSERT_NEAR(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_),10000*VTest::releps());
  }
  
  void do_spmv_test_single(double alpha, double beta, phist_d_complex sigma, phist_ZsparseMat_ptr z_A_shift)
  {
  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_r_[i]=ct::real(sigma);;
    x_A_->sigma_i_[i]=ct::imag(sigma);;
  }

  // sanity check of initial status
  ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
  ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_));

  phist_d_complex z_alpha=(phist_d_complex)alpha;
  phist_d_complex z_beta=(phist_d_complex)beta;

  SUBR(x_sparseMat_times_mvec)(alpha, x_A_, x_vec1_, beta, x_vec2_, &iflag_);
  ASSERT_EQ(0,iflag_);

  // check against multiplying with A-sigma*I directly
  phist_ZsparseMat_times_mvec(z_alpha,z_A_shift,z_vec1_,z_beta,z_vec2_,&iflag_);
  ASSERT_EQ(0,iflag_);

    ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
    ASSERT_NEAR(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_),1.0e-9);
  }
  
  // identity matrix (only used for checking if CARP is implemented at all right now)
  TYPE(sparseMat_ptr) I_;
  // backup vectors since the carp kernel works in-place
  TYPE(mvec_ptr) vec1b_,vec2b_,vec3b_;
  // mere pointers to allow e.g. passing in b=NULL
  TYPE(mvec_ptr) x_r, x_i, x_r_bak, x_i_bak, b;

  _MT_ sigma_r_[_NV_], sigma_i_[_NV_], omega_[_NV_];
  _CT_ sigma_[_NV_], minus_sigma_[_NV_];

  bool carpImplemented_, rc_carpImplemented_;
  bool cTypeImplemented_;
  
  // data structures for testing the "x_*" kernels used to implement CARP-CG.
  // these use x_mvecs with separate real and imaginary part so that CARP-CG can
  // be used in real arithmetic but with complex-shifted matrix
  TYPE(x_sparseMat) *x_A_ = NULL;
  TYPE(x_mvec) *x_vec1_ = NULL, *x_vec2_ = NULL,*x_vec3_ = NULL;
  
  phist_ZsparseMat_ptr z_A_, z_A_shift0_, z_A_shift1_;
  phist_Zmvec_ptr z_vec1_, z_vec2_,z_vec3_;
};

  TEST_F(CLASSNAME, create_matrices)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      ASSERT_TRUE(AssertNotNull(A_));
      if (cTypeImplemented_) 
      {
        ASSERT_TRUE(AssertNotNull(z_A_));
        ASSERT_TRUE(AssertNotNull(z_A_shift0_));
        ASSERT_TRUE(AssertNotNull(z_A_shift1_));
      }
      ASSERT_TRUE(AssertNotNull(I_));
    }
  }


  // make sure the operator is deterministic (important for CGMN)
  TEST_F(CLASSNAME, rc_operator_is_deterministic)
  {
    if (typeImplemented_ && !problemTooSmall_ && rc_carpImplemented_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      for (int i=0; i<_NV_; i++)
      {
        sigma_r_[i]=(MT)i;
        sigma_i_[i]=mt::one()-mt::one()/(MT)(i+1);
        omega_[i]=mt::one()+mt::one()/(MT)(i+1);
      }
      
      x_r=vec1_; x_r_bak=vec1b_;
      x_i=vec2_; x_i_bak=vec2b_;
      b=vec3_;
      create_and_apply_carp_rc(A_);
      ASSERT_EQ(0,iflag_);
            
      MT norm1r[_NV_], norm1i[_NV_];
      SUBR(mvec_norm2)(x_r,norm1r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_norm2)(x_i,norm1i,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // reset x_r and x_i to original vectors
      SUBR(mvec_add_mvec)(st::one(),vec1b_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec2b_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      create_and_apply_carp_rc(A_);
      ASSERT_EQ(0,iflag_);

      MT norm2r[_NV_], norm2i[_NV_];
      SUBR(mvec_norm2)(x_r,norm2r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_norm2)(x_i,norm2i,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MT_Test::ArraysEqual(norm1r,norm2r,_NV_,1,_NV_,1),10*mt::eps());
      ASSERT_NEAR(mt::one(),MT_Test::ArraysEqual(norm1i,norm2i,_NV_,1,_NV_,1),10*mt::eps());
    }
  }

  TEST_F(CLASSNAME, operator_is_deterministic)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      for (int i=0; i<_NV_; i++)
      {
        sigma_r_[i]=(MT)i;
      }
      
      x_r=vec1_; x_r_bak=vec1b_;
      b=vec3_;
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);
            
      MT norm1r[_NV_];
      SUBR(mvec_norm2)(x_r,norm1r,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // reset x_r and x_i to original vectors
      SUBR(mvec_add_mvec)(st::one(),vec1b_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);

      MT norm2r[_NV_];
      SUBR(mvec_norm2)(x_r,norm2r,&iflag_);

      ASSERT_NEAR(mt::one(),MT_Test::ArraysEqual(norm1r,norm2r,_NV_,1,_NV_,1),10*mt::eps());
    }
  }


  // test if the kernel works correctly if b=NULL is given (should be same as b=zeros(n,1))
  // and check if the operator is symmetric. Disabled for now because the symmetry-check seems
  // to fail with our implementations, I need to figure out if it is correct.
  TEST_F(CLASSNAME, DISABLED_operator_symmetric)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(vec2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(vec3_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      
      x_r=vec1_; x_r_bak=vec1b_;
      x_i=vec2_; x_i_bak=vec2b_;
      b=vec3_;
      for (int i=0; i<nvec_; i++)
      {
        sigma_i_[i]=mt::zero();
      }
      // copies x_r_bak=x_r, x_i_bak=x_i before the carp sweep
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);

      check_symmetry(x_r_bak,x_r,100*releps(x_r_bak));
      
      // check that x_i is still 0 (shift was s[j]+0i)
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec2_,st::zero()));
      
      x_r=vec2_; x_r_bak=vec2b_;
      x_i=NULL; x_i_bak=NULL;
      b=NULL;
      // reset x_r to original values, and try with b=x_i=NULL (should give same result)
      SUBR(mvec_add_mvec)(st::one(),vec1b_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec1_,st::one(),vec2_,&iflag_);
      ASSERT_NEAR(mt::one(),MvecEqual(vec2_,mt::zero()),10*VTest::releps());
    }
  }

  // test if the kernel works correctly if b=NULL is given (should be same as b=zeros(n,1))
  TEST_F(CLASSNAME, DISABLED_rc_operator_hermitian)
  {
    if (typeImplemented_ && !problemTooSmall_ && rc_carpImplemented_)
    {
      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(vec3_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      
      x_r=vec1_; x_r_bak=vec1b_;
      x_i=vec2_; x_i_bak=vec2b_;
      b=vec3_;
      // copies x_r_bak=x_r, x_i_bak=x_i before the carp sweep
      create_and_apply_carp_rc(A_);
      ASSERT_EQ(0,iflag_);
      
      check_symmetry_rc(x_r_bak,x_i_bak,x_r,x_i,10*releps(x_r_bak));
    }
  }

#if MATNAME==MATNAME_IDFUNC

  TEST_F(CLASSNAME, ZsparseMat_times_mvec_works)
  {
    if (!typeImplemented_ || problemTooSmall_ || !cTypeImplemented_) return;
    // the matrix z_A_shift0 is I - (1-i)I=I*i, so it's effect is to swap real and imagineary part and reverse the sign 
    // of the (new) real part.
    ASSERT_REAL_EQ(sigma_r_[0],1.0);
    ASSERT_REAL_EQ(sigma_i_[0],-1.0);
    phist_Zmvec_random(z_vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_Zmvec_add_mvec(ct::cmplx_I(),z_vec1_,ct::zero(),z_vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_ZsparseMat_times_mvec(ct::one(),z_A_shift0_,z_vec1_,-ct::one(),z_vec2_,&iflag_);
    MT nrm0[nvec_];
    phist_Zmvec_norm2(z_vec2_,nrm0,&iflag_);
    ASSERT_NEAR(mt::one(),ArrayEqual(nrm0,nvec_,1,nvec_,1,st::zero()),VTest::releps());
  }
  

  TEST_F(CLASSNAME, Identity_yields_zero)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      SUBR(mvec_put_value)(b,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(x_r,&iflag_);
      ASSERT_EQ(0,iflag_);
      for (int i=0; i<nvec_; i++)
      {
        sigma_r_[i]=st::zero();
        sigma_[i]=ct::zero();
        omega_[i]=1.0;
      }

      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);
      
      // real X should be zero now (all rows of I projected out)
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec1_,st::zero()));
      
    }
  }
#endif
//TODO - more tests like this could be invented, in particular
//       involving real or complex shifts

//////////////////////////////////////////////////
// tests for matrix/vector operations in "RC"   //
//////////////////////////////////////////////////

// these tests only use basic kernels (not the CARP kernel) but require complex
// arithmetic in order to work. So right now, the tests above will be run with
// builtin and epetra, and the tests below with ghost and tpetra only.
//
// At the start of these tests, z_vec1 = vec1 + i*vec1b etc

TEST_F(CLASSNAME, x_mvec_add_mvec)
{
  if (!cTypeImplemented_) return;

  double alpha = st::prand();
  double beta = st::prand();

  // sanity check of initial status
  ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
  ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_));
  
  phist_d_complex z_alpha = (phist_d_complex)alpha;
  phist_d_complex z_beta  = (phist_d_complex)beta;

  SUBR(x_mvec_add_mvec)(alpha,x_vec1_, beta, x_vec2_, &iflag_);
  ASSERT_EQ(0,iflag_);

  phist_Zmvec_add_mvec(z_alpha,z_vec1_, z_beta, z_vec2_, &iflag_);
  ASSERT_EQ(0,iflag_);
    
  ASSERT_NEAR(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_),mt::eps());
  ASSERT_NEAR(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_),1000*VTest::releps());
  
}

TEST_F(CLASSNAME, x_mvec_vadd_mvec)
{
  if (!cTypeImplemented_) return;

  double alpha[_NV_], alpha_i[_NV_];
  phist_d_complex z_alpha[_NV_];
  for (int i=0; i<nvec_;i++)
  {
    z_alpha[i]=ct::prand();
    alpha[i] = ct::real(z_alpha[i]);
    alpha_i[i] = ct::imag(z_alpha[i]);
  }
  double beta = st::prand();
  phist_d_complex z_beta  = (phist_d_complex)beta;

  // sanity check of initial status
  ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
  ASSERT_REAL_EQ(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_));
  

  SUBR(x_mvec_vadd_mvec)(alpha,alpha_i,x_vec1_, beta, x_vec2_, &iflag_);
  ASSERT_EQ(0,iflag_);

  phist_Zmvec_vadd_mvec(z_alpha,z_vec1_, z_beta, z_vec2_, &iflag_);
  ASSERT_EQ(0,iflag_);
    
  ASSERT_NEAR(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_),mt::eps());
  ASSERT_NEAR(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_),1.0e-9);
  
}

TEST_F(CLASSNAME, x_mvec_dot_mvec)
{
  if (!cTypeImplemented_) return;

  // sanity check of initial status
  ASSERT_EQ(1.0,MvecsEqualZD(z_vec1_, x_vec1_->v_, x_vec1_->vi_));
  ASSERT_EQ(1.0,MvecsEqualZD(z_vec2_, x_vec2_->v_, x_vec2_->vi_));
  
  double d_dot_r[nvec_], d_dot_i[nvec_];
  phist_d_complex z_dot[nvec_];

  SUBR(x_mvec_dot_mvec)(x_vec1_, x_vec2_, d_dot_r, d_dot_i, &iflag_);
  ASSERT_EQ(0,iflag_);

  phist_Zmvec_dot_mvec(z_vec1_, z_vec2_, z_dot, &iflag_);
  ASSERT_EQ(0,iflag_);

  double max_err_r=0.0, max_err_i=0.0;
  for (int i=0; i<nvec_; i++)
  {
    max_err_r = std::max(max_err_r,d_dot_r[i]-ct::real(z_dot[i]));
    max_err_i = std::max(max_err_i,d_dot_i[i]-ct::imag(z_dot[i]));
  }
  
  ASSERT_NEAR(1.0,1.0+max_err_r,std::sqrt(st::eps()));
  ASSERT_NEAR(1.0,1.0+max_err_i,std::sqrt(st::eps()));
  
}


TEST_F(CLASSNAME, x_sparseMat_times_mvec_without_shift)
{
  if (!cTypeImplemented_) return;

  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_r_[i]=0.0;
    x_A_->sigma_i_[i]=0.0;
  }
//  SUBR(mvec_put_value)(x_vec2_->v_,st::zero(),&iflag_);
//  SUBR(mvec_put_value)(x_vec2_->vi_,st::one(),&iflag_);
//  phist_Zmvec_put_value(z_vec2_,ct::cmplx_I(),&iflag_);
  do_spmv_test(1.0,0.0);
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_without_shift_with_alpha_beta)
{
  if (!cTypeImplemented_) return;

  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_r_[i]=0.0;
    x_A_->sigma_i_[i]=0.0;
  }
  do_spmv_test(1.234,-0.987);
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_real_shift)
{
  if (!cTypeImplemented_) return;

  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_r_[i]=0.0;
  }
  do_spmv_test(1.0,0.0);
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_real_shift_with_alpha_beta)
{
  if (!cTypeImplemented_) return;

  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_r_[i]=0.0;
  }
  do_spmv_test(1.234,-0.987);
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_imag_shift)
{
  if (!cTypeImplemented_) return;

  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_i_[i]=0.0;
  }
  do_spmv_test(1.0,0.0);
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_imag_shift_with_alpha_beta)
{
  if (!cTypeImplemented_) return;

  for (int i=0; i<nvec_; i++)
  {
    x_A_->sigma_i_[i]=0.0;
  }
  do_spmv_test(1.234,-0.987);
}


TEST_F(CLASSNAME, x_sparseMat_times_mvec_compare_with_i_shifted_matrix)
{
  if (!cTypeImplemented_) return;
  
  do_spmv_test_single(1.0,0.0, sigma_[0], z_A_shift0_);  
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_compare_with_rnd_shifted_matrix)
{
  if (!cTypeImplemented_) return;
  
  do_spmv_test_single(1.0,0.0, sigma_[1], z_A_shift1_);  
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_compare_with_i_shifted_matrix_alpha_beta)
{
  if (!cTypeImplemented_) return;
  
  double alpha = 1.23456;
  double beta = 9.876543;
  do_spmv_test_single(alpha, beta, sigma_[0], z_A_shift0_);  
}

TEST_F(CLASSNAME, x_sparseMat_times_mvec_compare_with_rnd_shifted_matrix_alpha_beta)
{
  if (!cTypeImplemented_) return;
  
  double alpha = 1.23456;
  double beta = 9.876543;
  do_spmv_test_single(alpha, beta, sigma_[1], z_A_shift1_);  

}

#endif
