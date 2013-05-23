#ifndef PHIST_KERNEL_TEST_WITH_VECTORS_H
#define PHIST_KERNEL_TEST_WITH_VECTORS_H

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "kernels/phist_kernels.h"
#include "KernelTestWithMap.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

using namespace testing;

/*! Test fixure. */
template<int _Nglob, int _Nvec>
class KernelTestWithVectors: public KernelTestWithMap<_Nglob> 
  {

public:

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithMap<_Nglob>::SetUp();
    nvec_=_Nvec;
    if (this->haveS_)
      {
      phist_Smvec_create(&svec1_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Smvec_extract_view(svec1_,&svec1_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Smvec_create(&svec2_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Smvec_extract_view(svec2_,&svec2_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      }
    if (this->haveD_)
      {
      phist_Dmvec_create(&dvec1_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Dmvec_extract_view(dvec1_,&dvec1_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);

      phist_Dmvec_create(&dvec2_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Dmvec_extract_view(dvec2_,&dvec2_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      }
    if (this->haveC_)
      {
      phist_Cmvec_create(&cvec1_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Cmvec_extract_view(cvec1_,&cvec1_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);

      phist_Cmvec_create(&cvec2_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Cmvec_extract_view(cvec2_,&cvec2_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      }
    if (this->haveZ_)
      {
      phist_Zmvec_create(&zvec1_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Zmvec_extract_view(zvec1_,&zvec1_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);

      phist_Zmvec_create(&zvec2_,this->map_,nvec_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      phist_Zmvec_extract_view(zvec2_,&zvec2_vp_,&lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithMap<_Nglob>::TearDown();
    if (this->haveS_)
      {
      phist_Smvec_delete(svec1_,&this->ierr_);
      phist_Smvec_delete(svec2_,&this->ierr_);
      }
    if (this->haveD_)
      {
      phist_Dmvec_delete(dvec1_,&this->ierr_);
      phist_Dmvec_delete(dvec2_,&this->ierr_);
      }
    if (this->haveC_)
      {
      phist_Cmvec_delete(cvec1_,&this->ierr_);
      phist_Cmvec_delete(cvec2_,&this->ierr_);
      }
    if (this->haveZ_)
      {
      phist_Zmvec_delete(zvec1_,&this->ierr_);
      phist_Zmvec_delete(zvec2_,&this->ierr_);
      }
    }

  int nloc_, nvec_, lda_;
  
    Smvec_ptr_t svec1_, svec2_;
    float *svec1_vp_, *svec2_vp_;
    
  Dmvec_ptr_t dvec1_, dvec2_;
  double *dvec1_vp_, *dvec2_vp_;

  Cmvec_ptr_t cvec1_, cvec2_;
  std::complex<float> *cvec1_vp_, *cvec2_vp_;

  Zmvec_ptr_t zvec1_, zvec2_;
  std::complex<double> *zvec1_vp_, *zvec2_vp_;


  };

#endif
