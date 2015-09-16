#include "../tools/TestHelpers.h"

#ifndef CLASSNAME
#error "file not included correctly."
#endif

extern "C" {

// prototype for a useful function from the driver_utils, we can't include
// the header here because it can only be included once after a phist_gen_X header.
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, const_comm_ptr_t comm,
        const char* problem, int* iflag);

} //extern "C"

/*! Test fixure. 
  
  basic tests for CARP kernel, creates a test matrix defined by
  _MATNAME_ and some vectors, checks if the kernel
  works with and without RHS, with and without real or complex shift.
  
  It is not easy to rigorously test this kernel routine because we give
  the kernel lib quite some freedom in deciding e.g. on the order of Kaczmarz
  updates and the parallelization scheme. What we can test is
  
  * every vector element should be updated
  * operator should be symmetric (hermitian for complex shift), S=X'(DKSWP(X))=S'
  
  given an Identity matrix,
  
  * for rhs 0 and shift 0  x=>0
  * for rhs 0 and shift 1  x=>x
  ... TODO: add tests ...
  
*/
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_> 
{

  public:

  typedef KernelTestWithType<MT> MT_Test;

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    // created in rebuildVectors
    vec1b_=NULL;
    vec2b_=NULL; 
    vec3b_=NULL;
    A_=NULL;
    I_=NULL;
    if (typeImplemented_ && !problemTooSmall_)
    {
      iflag_=PHIST_SPARSEMAT_OPT_CARP;
      SUBR(create_matrix)(&A_,comm_,_MATNAME_,&iflag_);
      ASSERT_EQ(0,iflag_);
      iflag_=PHIST_SPARSEMAT_OPT_CARP;
      SUBR(sparseMat_create_fromRowFunc)(&I_,comm_,_N_,_N_,1,NULL,&SUBR(idfunc),&iflag_);
      ASSERT_EQ(0,iflag_);

      for (int i=0; i<_NV_; i++)
      {
        sigma_r_[i]=mt::zero();
        sigma_i_[i]=mt::zero();
        omega_[i]=mt::one();
      }
    
      // check if CARP is implemented at all:
      void* work;
      SUBR(carp_setup)(I_, 1, sigma_r_, sigma_i_,
          &work, &iflag_);
      if (iflag_==PHIST_NOT_IMPLEMENTED) 
      {
        carpImplemented_=false;
      }
      else 
      {
        carpImplemented_=true;
        SUBR(carp_destroy)(I_, work, &iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    if (typeImplemented_ && !problemTooSmall_)
    {
      if (vec1b_!=NULL)
      {
        SUBR(mvec_delete)(vec1b_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if (vec2b_!=NULL)
      {
        SUBR(mvec_delete)(vec2b_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      if (vec3b_!=NULL)
      {
        SUBR(mvec_delete)(vec3b_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
      ASSERT_EQ(0,delete_mat(A_));
      ASSERT_EQ(0,delete_mat(I_));
    }
  }

// the matrices may have individual maps, so we need to recreate all vectors with the specific map of the matrix!
void rebuildVectors(TYPE(const_sparseMat_ptr) A)
{
  if (typeImplemented_ && !problemTooSmall_)
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

    if (vec1b_)
    {
      SUBR(mvec_delete)(vec1b_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    PHISTTEST_MVEC_CREATE(&vec1b_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    if (vec2b_)
    {
      SUBR(mvec_delete)(vec2b_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    PHISTTEST_MVEC_CREATE(&vec2b_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    if (vec3b_)
    {
      SUBR(mvec_delete)(vec3b_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    PHISTTEST_MVEC_CREATE(&vec3b_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);

    phist_map_get_local_length(domain_map, &nloc_, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    // set pointers for the tests
    x_r=vec1_; x_r_bak=vec1b_;
    x_i=vec2_; x_i_bak=vec2b_;
    b=vec3_;

  }
}


  // helper function to apply [x_r, x_i] = dkswp(A-sigma[j]*I, b, x_r, x_i)
  // and copy original X vectors to x_r_bak, x_i_bak.
  void create_and_apply_carp(TYPE(const_sparseMat_ptr) A)
  {

    SUBR(mvec_add_mvec)(st::one(),x_r,st::zero(),x_r_bak,&iflag_);
    ASSERT_EQ(0,iflag_);

    if (x_i!=NULL)
    {
      SUBR(mvec_add_mvec)(st::one(),x_i,st::zero(),x_i_bak,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    void* work;
    SUBR(carp_setup)(A, _NV_, sigma_r_, sigma_i_, &work, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    SUBR(carp_sweep)(A, sigma_r_, sigma_i_,b,x_r,x_i, work, omega_, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    SUBR(carp_destroy)(A, work, &iflag_);
    ASSERT_EQ(0,iflag_);
    
    return;
  }

protected:

int delete_mat(TYPE(sparseMat_ptr) &A)
  {
  if (A!=NULL)
    {
    SUBR(sparseMat_delete)(A,&iflag_);
    A = 0;
    }
  return iflag_;
  }

void check_symmetry(TYPE(const_mvec_ptr) X, TYPE(const_mvec_ptr) OPX,_MT_ tol=10*mt::eps())
  {
    _MT_ max_err=mt::zero();
    if (typeImplemented_ && !problemTooSmall_)
    {
      TYPE(sdMat_ptr) M=NULL;
      SUBR(sdMat_create)(&M,_NV_,_NV_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvecT_times_mvec)(st::one(),X,OPX,st::zero(),M,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_from_device)(M,&iflag_);
      ASSERT_EQ(0,iflag_);
      lidx_t ldm;
      _ST_* M_raw=NULL;
      SUBR(sdMat_extract_view)(M,&M_raw,&ldm,&iflag_);
      for (int i=0; i<_NV_; i++)
      {
        for (int j=i+1; j<_NV_;j++)
        {
          max_err = std::max(max_err,std::abs(M_raw[i*ldm+j]-M_raw[j*ldm+i]));
        }
      }
      SUBR(sdMat_delete)(M,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
    ASSERT_NEAR(mt::one(),max_err+mt::one(),tol);
  }

  TYPE(sparseMat_ptr) A_, I_;
  TYPE(mvec_ptr) vec1b_,vec2b_,vec3b_;
  // mere pointers to allow e.g. passing in b=NULL
  TYPE(mvec_ptr) x_r, x_i, x_r_bak, x_i_bak, b;

  _MT_ sigma_r_[_NV_], sigma_i_[_NV_], omega_[_NV_];

  bool carpImplemented_;

};

  TEST_F(CLASSNAME, create_matrices)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      ASSERT_TRUE(AssertNotNull(A_));
      ASSERT_TRUE(AssertNotNull(I_));
    }
  }


  // make sure the operator is deterministic (important for CGMN)
  TEST_F(CLASSNAME, operator_is_deterministic)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      rebuildVectors(A_);

      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      for (int i=0; i<_NV_; i++)
      {
        sigma_r_[i]=(MT)i;
        sigma_i_[i]=mt::one()-mt::one()/(MT)(i+1);
        omega_[i]=mt::one()+mt::one()/(MT)(i+1);
      }
      
      x_r=vec1_; x_r_bak=vec1b_;
      x_i=vec2_; x_i_bak=vec2b_;
      b=vec3_;
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);
            
      MT norm1r[_NV_], norm1i[_NV_];
      SUBR(mvec_norm2)(x_r,norm1r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_norm2)(x_i,norm1i,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // reset x_r and x_i to original vectors
      SUBR(mvec_add_mvec)(st::one(),vec1b_,st::zero(),vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec2b_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);

      MT norm2r[_NV_], norm2i[_NV_];
      SUBR(mvec_norm2)(x_r,norm2r,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_norm2)(x_i,norm2i,&iflag_);
      ASSERT_EQ(0,iflag_);

      ASSERT_NEAR(mt::one(),MT_Test::ArraysEqual(norm1r,norm2r,_NV_,1,_NV_,1),10*mt::eps());
      ASSERT_NEAR(mt::one(),MT_Test::ArraysEqual(norm1i,norm2i,_NV_,1,_NV_,1),10*mt::eps());
    }
  }


  // test if the kernel works correctly if b=NULL is given (should be same as b=zeros(n,1))
  TEST_F(CLASSNAME, works_with_bnull)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      rebuildVectors(A_);

      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(vec2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(vec3_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);
      
      x_r=vec1_; x_r_bak=vec1b_;
      x_i=vec2_; x_i_bak=vec2b_;
      b=vec3_;
      // copies x_r_bak=x_r, x_i_bak=x_i before the carp sweep
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);
      
      // check that x_i is still 0 (shift was 0+0i)
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec2_,st::zero()));
      
      x_r=vec2_; x_r_bak=vec2b_;
      x_i=NULL; x_i_bak=NULL;
      b=NULL;
      // reset x_r to original values, and try with b=x_i=NULL (should give same result)
      SUBR(mvec_add_mvec)(st::one(),vec1b_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      create_and_apply_carp(A_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(-st::one(),vec1_,st::one(),vec2_,&iflag_);
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec2_,mt::zero()));
      check_symmetry(x_r_bak,x_r,10*releps(x_r_bak));
    }
  }

  TEST_F(CLASSNAME, Identity_yields_zero)
  {
    if (typeImplemented_ && !problemTooSmall_ && carpImplemented_)
    {
      // matrices may have different maps
      rebuildVectors(I_);
      
      SUBR(mvec_random)(x_r,&iflag_);
      ASSERT_EQ(0,iflag_);

      create_and_apply_carp(I_);
      ASSERT_EQ(0,iflag_);
      
      // real X should be zero now (all rows of I projected out)
      ASSERT_REAL_EQ(mt::one(),MvecEqual(vec1_,st::zero()));
      
    }
  }

//TODO - more tests like this could be invented, in particular
//       involving real or complex shifts
