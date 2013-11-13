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
      read_mat(MATNAME,&A_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    if (typeImplemented_)
      {
      ASSERT_EQ(0,delete_mat(A_));
      }
    }

static const int numRuns_ = (int)(4.0e6/_N_/_NV_);
TYPE(crsMat_ptr) A_; 

protected:

int read_mat(const char* filebase,TYPE(crsMat_ptr) *ptr)
  {
  *ptr = NULL;
  char mmfile[256],hbfile[256],binfile[256];
  sprintf(mmfile,"%s.mm",filebase);
#ifdef IS_COMPLEX
  sprintf(hbfile,"%s.cua",filebase);
#else  
  sprintf(hbfile,"%s.rua",filebase);
#endif
  sprintf(binfile,"%s.bin",filebase);
  
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

};

  TEST_F(CLASSNAME, read_matrix) 
    {
    if (typeImplemented_)
      {
      ASSERT_TRUE(AssertNotNull(A_));
      }
    }

  TEST_F(CLASSNAME, times_mvec) 
    {
    if (typeImplemented_)
      {
      SUBR(mvec_random)(vec1_,&ierr_);
      SUBR(mvec_random)(vec2_,&ierr_);
      PHIST_SOUT(PHIST_INFO, "running crsMat_times_mvec %d times", numRuns_);
      for (int i=0; i<numRuns_;i++)
        {
        SUBR(crsMat_times_mvec)(1.0,A_,vec1_,0.0,vec2_,&ierr_);
        ASSERT_EQ(0,ierr_);
        }
      }
    }


