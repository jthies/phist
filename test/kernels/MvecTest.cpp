
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

#define _N_ 6
#define _NV_ 3

/*! Test fixure. */
class MvecTest: public KernelTestWithVectors<_N_,_NV_> {

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

  TEST_F(MvecTest, my_length) {
    int nloc;
    phist_Dmvec_my_length(dvec1_,&nloc,&ierr_);
    ASSERT_EQ(ierr_, 0);
    ASSERT_EQ(nloc, nglob_); // parallel testing not supported, yet.    
  }

  TEST_F(MvecTest, num_vectors) {
    int nvec;
    phist_Dmvec_num_vectors(dvec1_,&nvec,&ierr_);
    ASSERT_EQ(ierr_, 0);
    ASSERT_EQ(nvec, nvec_);
  }


  /*! Test the comm_get_rank function. */
  TEST_F(MvecTest, put_value) {
        phist_Dmvec_put_value(dvec1_,42.0,&ierr_);
	ASSERT_EQ(ierr_, 0);
	double sum=0.0;
	for (int i=0; i<nglob_;i++)
	  {
	  for (int j=0;j<nvec_;j++)
	    {
            sum+=abs(dvec1_vp_[j*lda_+i]-42.0);
            }
          }
        ASSERT_DOUBLE_EQ(sum,0.0);
}

