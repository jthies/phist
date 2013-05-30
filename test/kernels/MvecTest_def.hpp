#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;


#ifndef CLASSNAME
#error 'file not included correctly'
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
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    }

};

  TEST_F(CLASSNAME, my_length)
    {
    if (typeImplemented_)
      {
      int nloc;
      _SUBR_(mvec_my_length)(vec1_,&nloc,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(nglob_, nloc); // parallel testing not supported, yet.    
      }
    }

  TEST_F(CLASSNAME, num_vectors) 
    {
    if (typeImplemented_)
      {
      int nvec;
      _SUBR_(mvec_num_vectors)(vec1_,&nvec,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(nvec_, nvec);
      }
    }


  TEST_F(CLASSNAME, put_value) 
    {
    if (typeImplemented_)
      {
      _ST_ val = (_ST_)42.0 + (_ST_)3.0*_Complex_I;
      _SUBR_(mvec_put_value)(vec1_,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ((_MT_)1.0,ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,val));
      }
    }

  TEST_F(CLASSNAME, dot_mvec)
    {
    if (typeImplemented_)
      {
      for (int i=0;i<nloc_*nvec_;i++)
        {
        vec1_vp_[i]=random_number();
        vec2_vp_[i]=one()/_CONJ_(vec1_vp_[i]);
        }
      _ST_* dots = new _ST_[nvec_];
      _SUBR_(mvec_dot_mvec)(vec1_,vec2_,dots,&ierr_);
      ASSERT_EQ(0,ierr_);
            
      _ST_ val = one() * (_ST_)nloc_;
      ASSERT_REAL_EQ((_MT_)1.0,ArrayEqual(dots,nvec_,1,nvec_,1,val));
      delete [] dots;
      }
    
    }
    

