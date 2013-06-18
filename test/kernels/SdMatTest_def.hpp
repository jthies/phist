#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithSdMats.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;


#ifndef CLASSNAME
#error 'file not included correctly'
#endif
/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_> 
  {

public:


  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_>::SetUp();
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithSdMats<_ST_,_NROWS_,_NCOLS_>::TearDown();
    }

};

  TEST_F(CLASSNAME, put_value) 
    {
    if (typeImplemented_)
      {
      _ST_ val = (_ST_)42.0 + (_ST_)3.0*_Complex_I;
      _SUBR_(sdMat_put_value)(mat1_,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ((_MT_)1.0,ArrayEqual(mat1_vp_,nrows_,ncols_,m_lda_,1,val));
      }
    }
    

