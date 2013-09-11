#ifndef CLASSNAME
#error "file not included correctly."
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
      read_mat("spzero",&A0_);
      read_mat("speye",&A1_);
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
      ASSERT_EQ(0,delete_mat(A0_));
      ASSERT_EQ(0,delete_mat(A1_));
      ASSERT_EQ(0,delete_mat(A2_));
      ASSERT_EQ(0,delete_mat(A3_));
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
#ifdef _IS_COMPLEX_
  sprintf(hbfile,"%s%s%d.cua",_TPC_,filebase,nglob_);
#else  
  sprintf(hbfile,"%s%s%d.rua",_TPC_,filebase,nglob_);
#endif
  sprintf(binfile,"%s%s%d.bin",_TPC_,filebase,nglob_);
  
//  std::cout << "Looking for matrix \'"<<filebase<<"\'...\n";
//  std::cout << "... try \'"<<mmfile<<"\'\n";
  SUBR(crsMat_read_mm)(ptr,mmfile,&ierr_);
  if (ierr_!=_PHIST_SUCCESS_) // kernel lib can't read MatrixMarket format or file not found
    {
//    std::cout << "... try \'"<<hbfile<<"\'\n";
    SUBR(crsMat_read_hb)(ptr,hbfile,&ierr_);
    if (ierr_!=_PHIST_SUCCESS_) // kernel lib can't read Harwell-Boeing or file not found
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
      SUBR(mvec_put_value)(vec1_,val,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(1.0,A2_,vec1_,0.0,vec2_,&ierr_);
      if (ierr_) return (_MT_)ierr_;
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

  TEST_F(CLASSNAME, A0_times_mvec) 
    {
    if (typeImplemented_ && haveMats_)
      {
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(1.0,A0_,vec1_,0.0,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,0.0));
      }
    }


  TEST_F(CLASSNAME, A1_times_mvec)
    {
    if (typeImplemented_ && haveMats_)
      {
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      SUBR(crsMat_times_mvec)(1.0,A1_,vec1_,0.0,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
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

