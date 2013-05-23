
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 8
#define _NV_ 1

/*! Test fixure. */
class VecTest: public KernelTestWithVectors<_N_,_NV_> {

public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithVectors<_N_,_NV_>::SetUp();
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_N_,_NV_>::TearDown();
    }

};

  TEST_F(VecTest, my_length) {
    int nloc;
    phist_Dmvec_my_length(dvec1_,&nloc,&ierr_);
    ASSERT_EQ(ierr_, 0);
    ASSERT_EQ(nloc, nglob_); // parallel testing not supported, yet.    
  }


/*! Test the comm_get_rank function. */
TEST_F(VecTest, put_value) 
  {
  phist_Dmvec_put_value(dvec1_,42.0,&ierr_);
  ASSERT_EQ(ierr_, 0);
  double sum=0.0;
  for (int i=0; i<nglob_;i++)
    {
    sum+=dvec1_vp_[i]-42.0;
    }
  ASSERT_DOUBLE_EQ(sum,0.0);
  }

