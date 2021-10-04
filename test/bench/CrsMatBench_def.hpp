/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_,0,2> 
  {

public:
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,2>::VTest;


  /*! Set up routine.
   */
  virtual void SetUp()
    {
    VTest::SetUp();
    
    if (typeImplemented_ && !problemTooSmall_)
      {
      read_mat(MATNAME,&A_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    VTest::TearDown();
    if (typeImplemented_ && !problemTooSmall_)
      {
      ASSERT_EQ(0,delete_mat(A_));
      }
    }

static const int numRuns_ = (int)(4.0e6/_N_/_NV_);
TYPE(sparseMat_ptr) A_; 

protected:

int read_mat(const char* filebase,TYPE(sparseMat_ptr) *ptr)
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
  SUBR(sparseMat_read_mm)(ptr,mmfile,&iflag_);
  if (iflag_!=PHIST_SUCCESS) // kernel lib can't read MatrixMarket format or file not found
    {
//    std::cout << "... try \'"<<hbfile<<"\'\n";
    SUBR(sparseMat_read_hb)(ptr,hbfile,&iflag_);
    if (iflag_!=PHIST_SUCCESS) // kernel lib can't read Harwell-Boeing or file not found
      {
//      std::cout << "... try \'"<<binfile<<"\'\n";
      SUBR(sparseMat_read_bin)(ptr,binfile,&iflag_);
      }
    }
  return iflag_;
  }

int delete_mat(TYPE(sparseMat_ptr) A)
  {
  if (A!=NULL)
    {
    SUBR(sparseMat_delete)(A,&iflag_);
    }
  return iflag_;
  }

};

  TEST_F(CLASSNAME, read_matrix) 
    {
    if (typeImplemented_ && !problemTooSmall_)
      {
      ASSERT_TRUE(AssertNotNull(A_));
      }
    }

  TEST_F(CLASSNAME, times_mvec) 
    {
    if (typeImplemented_ && !problemTooSmall_)
      {
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      PHIST_SOUT(PHIST_INFO, "running sparseMat_times_mvec %d times", numRuns_);
      for (int i=0; i<numRuns_;i++)
        {
        SUBR(sparseMat_times_mvec)(1.0,A_,vec1_,0.0,vec2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        }
      }
    }


