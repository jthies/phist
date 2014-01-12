#ifndef CLASSNAME
#error "file not included correctly."
#endif

// the case A=0 does not work
// for ghost because the binCRS format
// doesn't handle it correctly, it seems
#ifdef PHIST_KERNEL_LIB_GHOST
#ifndef SKIP_ZERO_MAT
#define SKIP_ZERO_MAT
#endif
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
      SUBR(read_mat)("speye",nglob_,&A1_,&ierr_);
#ifndef SKIP_ZERO_MAT
      SUBR(read_mat)("spzero",nglob_,&A0_,&ierr_);
#else
      A0_=A1_;
#endif
      SUBR(read_mat)("sprandn",nglob_,&A2_,&ierr_);
      SUBR(read_mat)("sprandn_nodiag",nglob_,&A3_,&ierr_);
      SUBR(read_mat)("spshift",nglob_,&A4_,&ierr_);
      
      if (A0_==NULL || A1_==NULL || A2_==NULL || A3_==NULL || A4_==NULL)
        {
        haveMats_=false;
        }
      else
        {
        haveMats_=true;
        }
      }
    }

  /*! Clean up.
   */
  virtual void TearDown()
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    if (typeImplemented_)
      {
#ifndef SKIP_ZERO_MAT
      ASSERT_EQ(0,delete_mat(A0_));
#endif      
      ASSERT_EQ(0,delete_mat(A1_));
      ASSERT_EQ(0,delete_mat(A2_));
      ASSERT_EQ(0,delete_mat(A3_));
      ASSERT_EQ(0,delete_mat(A4_));
      }
    }

void rebuildVectors(TYPE(const_crsMat_ptr) A)
  {
  if (haveMats_ && typeImplemented_)
    {
    // set vec1 to be a valid X, vec2 a valid Y in Y=AX
    const_map_ptr_t range_map, domain_map;
    SUBR(crsMat_get_range_map)(A,&range_map,&ierr_);
    ASSERT_EQ(0,ierr_);
    SUBR(crsMat_get_domain_map)(A,&domain_map,&ierr_);
    ASSERT_EQ(0,ierr_);
    SUBR(mvec_delete)(vec1_,&ierr_);
    ASSERT_EQ(0,ierr_);
    SUBR(mvec_delete)(vec2_,&ierr_);
    ASSERT_EQ(0,ierr_);

    SUBR(mvec_create)(&vec1_,domain_map,this->nvec_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    int lda;
    SUBR(mvec_extract_view)(vec1_,&vec1_vp_,&lda,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    ASSERT_EQ(lda,this->lda_);
    SUBR(mvec_create)(&vec2_,range_map,this->nvec_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(mvec_extract_view)(vec2_,&vec2_vp_,&lda,&this->ierr_);
        ASSERT_EQ(0,this->ierr_);
    ASSERT_EQ(lda,this->lda_);
    }
  }

TYPE(crsMat_ptr) A0_; // all zero matrix
TYPE(crsMat_ptr) A1_; // identity matrix
TYPE(crsMat_ptr) A2_; // general sparse matrix with nonzero diagonal
TYPE(crsMat_ptr) A3_; // general sparse matrix with some zeros on the diagonal
TYPE(crsMat_ptr) A4_; // orthogonal sparse matrix that "shifts" each value to the next row

protected:

int delete_mat(TYPE(crsMat_ptr) A)
  {
  if (A!=NULL)
    {
    SUBR(crsMat_delete)(A,&ierr_);
    }
  return ierr_;
  }

_MT_ const_row_sum_test(TYPE(crsMat_ptr) A)
  {
    if (typeImplemented_ && haveMats_)
      {
      _ST_ val = random_number();
      global_sum(&val,1,mpi_comm_);
      // build vectors with correct range- and domain map
      rebuildVectors(A);
      SUBR(mvec_put_value)(vec1_,val,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(st::one(),A2_,vec1_,st::zero(),vec2_,&ierr_);
      if (ierr_) return (_MT_)ierr_;
      SUBR(mvec_print)(vec1_,&ierr_);
      SUBR(mvec_print)(vec2_,&ierr_);
      return ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,val);
      }
  return mt::one();
  }

  bool haveMats_;
};

  TEST_F(CLASSNAME, read_matrices) 
    {
    if (typeImplemented_)
      {
      ASSERT_TRUE(AssertNotNull(A0_));
      ASSERT_TRUE(AssertNotNull(A1_));
      ASSERT_TRUE(AssertNotNull(A2_));
      ASSERT_TRUE(AssertNotNull(A3_));
      }
    }

#ifndef SKIP_ZERO_MAT
  TEST_F(CLASSNAME, A0_times_mvec) 
    {
    if (typeImplemented_ && haveMats_)
      {
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(st::one(),A0_,vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,0.0));
      }
    }
#endif

  TEST_F(CLASSNAME, A1_times_mvec)
    {
    if (typeImplemented_ && haveMats_)
      {
      ST alpha, beta;
      //I*X=X?
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(st::one(),A1_,vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif

      //alpha*I*X=alpha*X?
      alpha = random_number();
      beta=st::zero();
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_scale)(vec1_,alpha,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif

      //0*I*X+beta*Y = beta*Y? 
      alpha=st::zero(); 
      beta=random_number();
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
      SUBR(mvec_random)(vec1_,&ierr_); 
      SUBR(mvec_random)(vec2_,&ierr_); 
#if PHIST_OUTLEV>=PHIST_DEBUG
      std::cout << "input="<<std::endl;
      SUBR(mvec_print)(vec1_,&ierr_);
      std::cout << "output, before="<<std::endl;
      SUBR(mvec_print)(vec2_,&ierr_);
#endif
      // v3=beta*v2 
      SUBR(mvec_add_mvec)(beta,vec2_,st::zero(),vec3_,&ierr_); 
      ASSERT_EQ(0,ierr_); 
      // v2 = 0*v1 + beta*v2 (=v3) 
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_); 
      ASSERT_EQ(0,ierr_); 
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&ierr_);
      SUBR(mvec_print)(vec3_,&ierr_);
#endif
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_));
#endif

      //I*X+beta*Y = X+beta*Y?
      alpha = st::one();
      beta = random_number();
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      std::cout << "input="<<std::endl;
      SUBR(mvec_print)(vec1_,&ierr_);
      std::cout << "output, before="<<std::endl;
      SUBR(mvec_print)(vec2_,&ierr_);
#endif
      //v3=v1+beta*v2
      SUBR(mvec_add_mvec)(alpha,vec1_,st::zero(),vec3_,&ierr_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&ierr_);
      // v2 = v1 + beta*v2 (=alpha*v1+v3)
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&ierr_);
      SUBR(mvec_print)(vec3_,&ierr_);
#endif
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_));
#endif

      //alpha*I*X+beta*Y = alpha*X+beta*Y?
      alpha = random_number();
      beta = random_number();
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      std::cout << "input="<<std::endl;
      SUBR(mvec_print)(vec1_,&ierr_);
      std::cout << "output, before="<<std::endl;
      SUBR(mvec_print)(vec2_,&ierr_);
#endif
       // v3=alpha*v1+beta*v2
      SUBR(mvec_add_mvec)(alpha, vec1_,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // v2 = alpha*v1 + beta*v2 (=alpha*v1+v3)
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&ierr_);
      SUBR(mvec_print)(vec3_,&ierr_);
#endif
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_));
#endif
      }
    }

  TEST_F(CLASSNAME, A2_times_mvec)
    {
    // we allow a tolerance here because the matrices may have errors in the
    // last digit and we can't get the test to pass otherwise.
    ASSERT_NEAR(mt::one(),const_row_sum_test(A2_),100*mt::eps());
    }

  TEST_F(CLASSNAME, A3_times_mvec)
    {
    // we allow a tolerance here because the matrices may have errors in the
    // last digit and we can't get the test to pass otherwise.
    ASSERT_NEAR(mt::one(),const_row_sum_test(A3_),100*mt::eps());
    }

  TEST_F(CLASSNAME, shift_mvec)
  {
    if( typeImplemented_ )
    {
      _ST_ alpha = st::one();
      _ST_ beta = st::zero();
      int ilower = 0;
      phist_map_get_ilower(map_,&ilower,&ierr_);
      ASSERT_EQ(0,ierr_);

      // setup recognizable input
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec1_vp_[i*lda_+j] = (_ST_)(ilower+i + j*nglob_);
#else
          vec1_vp_[j*lda_+i] = (_ST_)(ilower+i + j*nglob_);
#endif
      }

      // apply our shift matrix
      std::cout << "MVM with A='shift', alpha=1, beta=0"<<std::endl;
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#endif
      SUBR(crsMat_times_mvec)(alpha,A4_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#endif

      // check result
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
#ifdef PHIST_KERNEL_LIB_FORTRAN
          ASSERT_REAL_EQ((ilower+i+1)%nglob_ + j*nglob_, st::real(vec2_vp_[i*lda_+j]));
#else
          ASSERT_REAL_EQ((ilower+i+1)%nglob_ + j*nglob_, st::real(vec2_vp_[j*lda_+i]));
#endif
      }
    }
  }



  TEST_F(CLASSNAME, A2_precalc_result)
  {
    if( typeImplemented_ )
    {
      _ST_ alpha = st::one();
      _ST_ beta = st::zero();
      int ilower = 0;
      phist_map_get_ilower(map_,&ilower,&ierr_);
      ASSERT_EQ(0,ierr_);

      // setup recognizable input
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec1_vp_[i*lda_+j] = (_ST_)(ilower+i + j*nglob_);
#else
          vec1_vp_[j*lda_+i] = (_ST_)(ilower+i + j*nglob_);
#endif
      }

      // apply our shift matrix
      std::cout << "MVM with A='rand', alpha=1, beta=0"<<std::endl;
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
#endif
      SUBR(crsMat_times_mvec)(alpha,A2_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
#endif
 
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
#else
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
#endif
      // check result
      for(int i = 0; i < nloc_; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
#ifdef PHIST_KERNEL_LIB_FORTRAN
          ASSERT_NEAR(st::real(precalc_result[(ilower+i)*nvec_+j]), st::real(vec2_vp_[i*lda_+j]), mt::sqrt(mt::eps()));
          ASSERT_NEAR(st::imag(precalc_result[(ilower+i)*nvec_+j]), st::imag(vec2_vp_[i*lda_+j]), mt::sqrt(mt::eps()));
#else
          ASSERT_NEAR(st::real(precalc_result[(ilower+i)*nvec_+j]), st::real(vec2_vp_[j*lda_+i]), mt::sqrt(mt::eps()));
          ASSERT_NEAR(st::imag(precalc_result[(ilower+i)*nvec_+j]), st::imag(vec2_vp_[j*lda_+i]), mt::sqrt(mt::eps()));
#endif
        }
      }
    }
  }

