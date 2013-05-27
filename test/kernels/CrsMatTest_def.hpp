#include "phist_macros.h"

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace ::testing;

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
      ASSERT_EQ(read_mat("spzero",&A0_),0);
      ASSERT_EQ(read_mat("speye",&A1_),0);
      ASSERT_EQ(read_mat("sprandn",&A2_),0);
      ASSERT_EQ(read_mat("sprandn_no_diag",&A3_),0);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    if (typeImplemented_)
      {
      ASSERT_EQ(delete_mat(A0_),0);
      ASSERT_EQ(delete_mat(A1_),0);
      ASSERT_EQ(delete_mat(A2_),0);
      ASSERT_EQ(delete_mat(A3_),0);
      }
    }

_TYPE_(crsMat_ptr) A0_; // all zero matrix
_TYPE_(crsMat_ptr) A1_; // identity matrix
_TYPE_(crsMat_ptr) A2_; // general sparse matrix with nonzero diagonal
_TYPE_(crsMat_ptr) A3_; // general sparse matrix with some zeros on the diagonal

protected:

int read_mat(const char* filebase,_TYPE_(crsMat_ptr) *ptr)
  {
  char mmfile[256],hbfile[256],binfile[256];
  sprintf(mmfile,"%s%s%d.mm",_TPC_,filebase,nglob_);
  sprintf(hbfile,"%s%s%d.rua",_TPC_,filebase,nglob_);
  sprintf(binfile,"%s%s%d.bin",_TPC_,filebase,nglob_);
  
  std::cout << "Looking for matrix \'"<<filebase<<"\'...\n";
  std::cout << "... try \'"<<mmfile<<"\'\n";
  _SUBR_(crsMat_read_mm)(ptr,mmfile,&ierr_);
  if (ierr_==-99 || ierr_==-1) // kernel lib can't read MatrixMarket format or file not found
    {
    std::cout << "... try \'"<<hbfile<<"\'\n";
    _SUBR_(crsMat_read_hb)(ptr,hbfile,&ierr_);
    if (ierr_==-99 || ierr_==-1) // kernel lib can't read Harwell-Boeing or file not found
      {
      std::cout << "... try \'"<<binfile<<"\'\n";
      _SUBR_(crsMat_read_bin)(ptr,binfile,&ierr_);
      }
    }
  return ierr_;
  }

int delete_mat(_TYPE_(crsMat_ptr) A)
  {
  if (A!=NULL)
    {
    _SUBR_(crsMat_delete)(A,&ierr_);
    }
  return ierr_;
  }

::testing::AssertionResult const_row_sum_test(_TYPE_(crsMat_ptr) A)
  {
    if (typeImplemented_)
      {
      _ST_ val = random_number();
      _SUBR_(mvec_put_value)(vec1_,val,&ierr_);
      _SUBR_(mvec_random)(vec2_,&ierr_);
      _SUBR_(crsMat_times_mvec)(1.0,A2_,vec1_,0.0,vec2_,&ierr_);
      if (ierr_) return AssertionFailure();
      return ArrayEqual(vec2_vp_,nloc_,val);
      }
  return AssertionSuccess();
  }
};

  TEST_F(CLASSNAME, A0_times_mvec) 
    {
    if (typeImplemented_)
      {
      _SUBR_(mvec_random)(vec1_,&ierr_);
      _SUBR_(mvec_random)(vec2_,&ierr_);
      _SUBR_(crsMat_times_mvec)(1.0,A0_,vec1_,0.0,vec2_,&ierr_);
      ASSERT_EQ(ierr_,0);
      ASSERT_TRUE(ArrayEqual(vec2_vp_,nloc_,0.0));
      }
    }


  TEST_F(CLASSNAME, A1_times_mvec)
    {
    if (typeImplemented_)
      {
      _SUBR_(mvec_random)(vec1_,&ierr_);
      _SUBR_(mvec_random)(vec2_,&ierr_);
      _SUBR_(crsMat_times_mvec)(1.0,A0_,vec1_,0.0,vec2_,&ierr_);
      ASSERT_EQ(ierr_,0);
      ASSERT_TRUE(ArraysEqual(vec1_vp_,vec2_vp_,nloc_));
      }
    }

  TEST_F(CLASSNAME, A2_times_mvec)
    {
    ASSERT_TRUE(const_row_sum_test(A2_));
    }

  TEST_F(CLASSNAME, A3_times_mvec)
    {
    ASSERT_TRUE(const_row_sum_test(A3_));
    }

