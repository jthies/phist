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
      read_mat("speye",&A1_);
#ifndef SKIP_ZERO_MAT
      read_mat("spzero",&A0_);
#else
      A0_=A1_;
#endif
      read_mat("sprandn",&A2_);
      read_mat("sprandn_nodiag",&A3_);
      
      if (A0_==NULL || A1_==NULL || A2_==NULL || A3_==NULL)
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

protected:

int read_mat(const char* filebase,TYPE(crsMat_ptr) *ptr)
  {
  *ptr = NULL;
  char mmfile[256],hbfile[256],binfile[256];
  sprintf(mmfile,"%s%s%d.mm",_TPC_,filebase,nglob_);
#ifdef IS_COMPLEX
  sprintf(hbfile,"%s%s%d.cua",_TPC_,filebase,nglob_);
#else  
  sprintf(hbfile,"%s%s%d.rua",_TPC_,filebase,nglob_);
#endif
  sprintf(binfile,"%s%s%d.bin",_TPC_,filebase,nglob_);
  
//  std::cout << "Looking for matrix \'"<<filebase<<"\'...\n";
//  std::cout << "... try \'"<<mmfile<<"\'\n";
  SUBR(crsMat_read_mm)(ptr,mmfile,&ierr_);
  if (ierr_!=PHIST_SUCCESS) // kernel lib can't read MatrixMarket format or file not found
    {
//    std::cout << "... try \'"<<hbfile<<"\'\n";
    SUBR(crsMat_read_hb)(ptr,hbfile,&ierr_);
    if (ierr_!=PHIST_SUCCESS) // kernel lib can't read Harwell-Boeing or file not found
      {
//      std::cout << "... try \'"<<binfile<<"\'\n";
      SUBR(crsMat_read_bin)(ptr,binfile,&ierr_);
      }
    }
  return ierr_;
  }

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
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));

      //alpha*I*X=alpha*X?
      alpha = random_number();
      beta=st::zero();
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_scale)(vec1_,alpha,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));

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
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_));

      //I*X+beta*Y = X+beta*Y?
      alpha = st::one();
      beta = random_number();
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      //v3=v1+beta*v2
      SUBR(mvec_add_mvec)(alpha,vec1_,st::zero(),vec3_,&ierr_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&ierr_);
      // v2 = v1 + beta*v2 (=alpha*v1+v3)
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_));

      //alpha*I*X+beta*Y = alpha*X+beta*Y?
      alpha = random_number();
      beta = random_number();
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      // v3=alpha*v1+beta*v2
      SUBR(mvec_add_mvec)(alpha, vec1_,st::zero(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&ierr_);
      ASSERT_EQ(0,ierr_);
      // v2 = alpha*v1 + beta*v2 (=alpha*v1+v3)
      SUBR(crsMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_));
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

