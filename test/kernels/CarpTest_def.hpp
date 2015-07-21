#include "../tools/TestHelpers.h"
/* we don't have CARP tests yet */
#if 0
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. 

\todo this test class has a lot in common with
CrsMatTest, we should have a common base class like
KernelTestWithCrsMatrices<_N_> that reads the matrices etc.
*/
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_> 
{

  public:

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    vec2bak_=NULL; // created in rebuildVectors
    vec3bak_=NULL;
    if (typeImplemented_)
    {
      SUBR(read_mat)("speye",nglob_,&A1_,&iflag_);
      SUBR(read_mat)("sprandn",nglob_,&A2_,&iflag_);
      SUBR(read_mat)("sprandn_nodiag",nglob_,&A3_,&iflag_);
      
      if (A1_==NULL || A2_==NULL || A3_==NULL || A4_==NULL)
      {
        haveMats_=false;
      }
      else
      {
        haveMats_=true;
      }
    }
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    if (typeImplemented_)
    {
      if (vec2bak_!=NULL)
      {
        SUBR(mvec_delete)(vec2bak_,&ierr_);
      }
      if (vec3bak_!=NULL)
      {
        SUBR(mvec_delete)(vec3bak_,&ierr_);
      }
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(0,delete_mat(A1_));
      ASSERT_EQ(0,delete_mat(A2_));
      ASSERT_EQ(0,delete_mat(A3_));
      PHIST
    }
  }

// the matrices may have individual maps, so we need to recreate all vectors with the specific map of the matrix!
void rebuildVectors(TYPE(const_sparseMat_ptr) A)
{
  if (typeImplemented_ && haveMats_)
  {
    // set vec1 to be a valid X, vec2 and vec3 a valid Y in Y=AX
    const_map_ptr_t range_map, domain_map;
    SUBR(sparseMat_get_range_map)(A,&range_map,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sparseMat_get_domain_map)(A,&domain_map,&iflag_);

    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);

    lidx_t lda;
    PHISTTEST_MVEC_CREATE(&vec1_,domain_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec1_,&vec1_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    PHISTTEST_MVEC_CREATE(&vec2_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec2_,&vec2_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    PHISTTEST_MVEC_CREATE(&vec3_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec3_,&vec3_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    PHISTTEST_MVEC_CREATE(&vec2bak_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    PHISTTEST_MVEC_CREATE(&vec3bak_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);

    phist_map_get_local_length(domain_map, &nloc_, &iflag_);
    ASSERT_EQ(0,iflag_);
  }
}


  // helper function to apply [X_r, X_i] = dkswp(A-sigma*I, B, X_r, X_i)
  // for hsift sigma_r+i*sigma_i and B in vec1_, X_r in vec2_ and X_i in vec3_.
  void create_and_apply_carp
        (TYPE(const_sparseMat_ptr) A, _MT_ sigma_r, _MT_ sigma_i, _MT_ omega, int* iflag)
  {
    if( !typeImplemented_ )
      return;

    _MT_* nrms_ai2i;
    void* work;
    SUBR(carp_setup)(A, 1, &sigma_r, sigma_i,
        &nrms_ai2i_, &work, iflag);
    ASSERT_EQ(0,*iflag);
    SUBR(carp_sweep)(A, 1,&sigma_r,&sigma_i,vec1_,&vec2_,&vec3_,
        nrm_ai2i, work, &omega, iflag);
    ASSERT_EQ(0,*iflag);
    SUBR(carp_destroy)(A, 1,
        nrms_ai2i, work, iflag)
    ASSERT_EQ(0,*iflag);
  }

TYPE(sparseMat_ptr) A1_; // identity matrix
TYPE(sparseMat_ptr) A2_; // general sparse matrix with nonzero diagonal
TYPE(sparseMat_ptr) A3_; // general sparse matrix with some zeros on the diagonal
TYPE(sparseMat_ptr) A4_; // orthogonal sparse matrix that "shifts" each value to the next row

protected:

int delete_mat(TYPE(sparseMat_ptr) A)
  {
  if (A!=NULL)
    {
    SUBR(sparseMat_delete)(A,&iflag_);
    }
  return iflag_;
  }

_MT_ check_symmetry(TYPE(sparseMat_ptr) A)
{
    if (typeImplemented_ && haveMats_)
    {
      //TODO test not implemented
    }
  return mt::one();
}

 // assuming the "old" X is in vec2bak_ + i*vec3bal_, and the
 // updated one after the sweep in vec2_+i*vec3_, check some 
 // basic properties.
 
 // norm of x less or equal after the sweep to before (note that we
 // project out some components!)
_MT_ check_nonincreasing(TYPE(sparseMat_ptr) A)
{
    if (typeImplemented_ && haveMats_)
    {
      _MT_ nrms_old[nvec_], nrms_new[nvec_];
      _MT_ imTim_old[nvec_], imTim_new[nvec_];

      SUBR(mvec_norm2)(vec2bak_,nrms_old,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_norm2)(vec3bak_,imTim_old,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_norm2)(vec2_,nrms_new,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_norm2)(vec3_,imTim_new,&ierr_);
      ASSERT_EQ(0,ierr_);

    for (int i=0;i<nvec_;i++)
    {  
      nrms_old[i]=mt::sqrt(nrms_old[i]*nrms_old[i] + imTim_old[i]*imTim_old[i])
      nrms_new[i]=mt::sqrt(nrms_new[i]*nrms_new[i] + imTim_new[i]*imTim_new[i])
      ASSERT_TRUE(nrms_r_new[i]<=nrms_old[i]);
    }
  }
  return mt::one();
}

  bool haveMats_;
};

  TEST_F(CLASSNAME, read_matrices)
  {
    if (typeImplemented_)
    {
      ASSERT_TRUE(AssertNotNull(A0_));
      ASSERT_TRUE(AssertNotNull(A1_));
      ASSERT_TRUE(AssertNotNull(A2_));
      ASSERT_TRUE(AssertNotNull(A3_));
    }
  }

  TEST_F(CLASSNAME, A1_unshifted)
  {
    if (typeImplemented_ && haveMats_)
    {
      // matrices may have different maps
      rebuildVectors(A1_);

      MT_ sigma_r=mt::zero(), sigma_i=mt::zero();

      //CARP(A=I,b=0,x)=x?

      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_put_value)(vec2_,st::zero(),&iflag_);
      SUBR(mvec_put_value)(vec3_,st::zero(),&iflag_);

      create_and_apply_carp(A1_, mt::zero(), mt::zero(), mt::one(), &iflag_);
      ASSERT_EQ(0,iflag_);
      
      ASSERT_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_));
      ASSERT_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_));
    }
  }







  TEST_F(CLASSNAME, A2_unshifted)
  {
    if (typeImplemented_ && haveMats_)
    {
      // matrices may have different maps
      rebuildVectors(A2_);

      MT_ sigma_r=mt::zero(), sigma_i=mt::zero();

      //CARP(A=I,b=0,x)=x?

      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(mvec_put_value)(vec3_,st::zero(),&iflag_);

      create_and_apply_carp(A2, mt::zero(), mt::zero(), mt::one(), &iflag_);
      ASSERT_EQ(0,iflag_);
      
      // real X should stay real
      ASSERT_EQ(mt::one(),ArrayEqual(vec3_vp_,nloc_,nvec_,lda_,stride_,st::zero(),vflag_));
      
      // check basic invariants
      ASSERT_EQ(mt::one(),check_symmetry());
      ASSERT_EQ(mt::one(),check_nonincreasing());
  }
}
#endif
