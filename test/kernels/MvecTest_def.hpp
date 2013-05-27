#include "phist_macros.h"

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
      ASSERT_EQ(ierr_, 0);
      ASSERT_EQ(nloc, nglob_); // parallel testing not supported, yet.    
      }
    }

  TEST_F(CLASSNAME, num_vectors) 
    {
    if (typeImplemented_)
      {
      int nvec;
      _SUBR_(mvec_num_vectors)(vec1_,&nvec,&ierr_);
      ASSERT_EQ(ierr_, 0);
      ASSERT_EQ(nvec, nvec_);
      }
    }


  TEST_F(CLASSNAME, put_value) 
    {
    if (typeImplemented_)
      {
      _ST_ val = (_ST_)42.0 + (_ST_)3.0*_Complex_I;
      _SUBR_(mvec_put_value)(vec1_,val,&ierr_);
      ASSERT_EQ(ierr_, 0);
      ASSERT_TRUE(ArrayEqual(vec1_vp_,nloc_,val));
      }
    }

  TEST_F(CLASSNAME, dot_mvec)
    {
    if (typeImplemented_)
      {
      for (int i=0;i<nloc_*nvec_;i++)
        {
        vec1_vp_[i]=random_number();
        vec2_vp_[i]=one()/vec1_vp_[i];
        }
      _ST_* dots = new _ST_[nvec_];
      _SUBR_(mvec_dot_mvec)(vec1_,vec2_,dots,&ierr_);
      ASSERT_EQ(ierr_, 0);
      
      std::cout << "I="<<_Complex_I<<std::endl;
      std::cout << std::endl;
      for (int i=0;i<nloc_;i++)
        {
        for (int j=0;j<nvec_;j++)
          {
          std::cout << vec1_vp_[lda_*j+i] << " ";
          }
        std::cout << std::endl;
        }
      std::cout << std::endl;
      for (int i=0;i<nloc_;i++)
        {
        for (int j=0;j<nvec_;j++)
          {
          std::cout << vec2_vp_[lda_*j+i] << " ";
          }
        std::cout << std::endl;
        }

      _ST_ val = one() * (_ST_)nloc_;
      
      std::cout << "expected: "<<val<<std::endl;
      std::cout << std::endl;
      for (int i=0;i<nvec_;i++) 
        {
        std::cout << dots[i]<<std::endl;
        }
      
      ASSERT_TRUE(ArrayEqual(dots,nvec_,val));
      delete [] dots;
      }
    
    }
    

