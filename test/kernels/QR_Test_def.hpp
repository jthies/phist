#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithVectors.h"
#include "KernelTestWithSdMats.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;


#ifndef CLASSNAME
#error 'file not included correctly'
#endif
/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_> 
  {

public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::SetUp();
    if (typeImplemented_)
      {
      _SUBR_(mvec_random)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<stride_*nrows_;i+=stride_)
          {
          vec2_vp_[j*lda_+i] = vec1_vp_[j*lda_+i];
          }
        }
      _SUBR_(mvec_QR)(vec2_,mat1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    KernelTestWithSdMats<_ST_,_NV_,_NV_>::TearDown();
    }

};

  // check if vectors are normalized correctly after QR factorization
  TEST_F(CLASSNAME, normality) 
    {
    if (typeImplemented_)
      {
      // see if all columns in vec2 have 2-norm 1
      _ST_ *norms = new _ST_[ncols_];
      for (int j=0;j<ncols_;j++)
        {
        _ST_ sum=zero();
        for (int i=0;i<stride_*nrows_;i+=stride_)
          {
          _ST_ val=vec2_vp_[j*lda_+i];
          sum+=val*_CONJ_(val); 
          }
        norms[j]=std::sqrt(sum);
        }
      ASSERT_REAL_EQ((_MT_)1.0,ArrayEqual(norms,nvec_,1,nvec_,1,one()));
      delete [] norms;
      }
    }


  // check if vectors are mutually orthogonal after QR factorization
  TEST_F(CLASSNAME, ortho) 
    {
    if (typeImplemented_)
      {
      int nsums=(nrows_*ncols_-nrows_)/2;
      _ST_ sums[nsums];
      int k=0;
      for (int j1=0;j1<ncols_;j1++)
      for (int j2=j1+1;j2<ncols_;j2++)
        {
        _ST_ sum=zero();
        for (int i=0;i<stride_*nrows_;i+=stride_)
          {
          _ST_ val1=vec2_vp_[j1*lda_+i];
          _ST_ val2=vec2_vp_[j2*lda_+i];
          sum+=val1*_CONJ_(val2); 
          }
        sums[k++]=sum;
        }
      ASSERT_REAL_EQ((_MT_)1.0,ArrayEqual(sums,nsums,1,nvec_,1,zero()));
      }
    }

