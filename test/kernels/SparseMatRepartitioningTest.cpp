#include "phist_kernels.h"
#include "../tools/MatrixIO.h"

#ifdef PHIST_HAVE_ESSEX_PHYSICS
#include "essex-physics/matfuncs.h"
#else
#include "matfuncs.h"
#endif

#include "phist_gen_d.h"
#include "phist_driver_utils.h"
#include "../tools/TestHelpers.h"
#include "TestWithType.h"

#include "KernelTestWithMap.h"

/* need to include this again because TestWithType defines all variants
   and thus invalidates phist_gen_d above
 */
#include "phist_gen_d.h"

#include <cstdio>

//#define _L_ 12
//#define _N_ 924

//#define _L_ 10
//#define _N_ 252

#define _L_ 8
#define _N_ 70

using namespace phist::testing;

/* This test constructs the spinSZ[L] matrix with and without the flag    */
/* PHIST_SPARSEMAT_PERM_GLOBAL, constructs a vector X with characteristic */
/* entries, and computes the spMVM in three different ways, comparing the */
/* results:                                                               */
/* - by using the row function and vector definition directly (does not   */
/*   require any communication or kernels)                                */
/* - with the standard map (linearly partitioned), does not require       */
/*   mvec_to_mvec to work                                                 */
/* - with the possibly redistributed matrix and mvec_to_mvec to get the   */
/*   result in the original map                                           */
/*                                                                        */
/* We use the row map of the matrices for all vectors in this test, as is */
/* commonly done in PHIST. This also tests that the kernel library sticks */
/* to the unwritten rule that row, range and domain map should typically  */
/* be the same for our applications. Having a test here is useful for     */
/* verifying this behavior, too.                                          */
class DSparseMatRepartTest: public virtual TestWithType<double>,
                            public KernelTestWithMap<_N_>
{

  public:
  
  // construct the matrices and vectors,
  // and initialize input and exact result vector for spMVM
  virtual void SetUp()
  {
    // construct the spinSZ[L] matrix with the default map
    char spin_label[9];
#if _L_>8
    sprintf(spin_label,"spinSZ%2d",_L_);
#else
    sprintf(spin_label,"spinSZ%1d",_L_);
#endif
    iflag_=0;
    phist_Dcreate_matrix(&linearA_,comm_,spin_label,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_DsparseMat_get_domain_map(linearA_,&linearMap_,&iflag_);
    EXPECT_EQ(0,iflag_);
    // construct same matrix but possibly with repartitioning (depends on what the kernel
    // library supports and which TPLs are available)
    iflag_=PHIST_SPARSEMAT_PERM_GLOBAL|getSparseMatCreateFlag(_N_,1);
    phist_Dcreate_matrix(&repartA_,comm_,spin_label,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_DsparseMat_get_domain_map(repartA_,&repartMap_,&iflag_);
    EXPECT_EQ(0,iflag_);
    
    // create the vectors
    PHISTTEST_MVEC_CREATE(&linearV_in_,linearMap_,1,&iflag_);
    EXPECT_EQ(0,iflag_);
    PHISTTEST_MVEC_CREATE(&linearV_out_,linearMap_,1,&iflag_);
    EXPECT_EQ(0,iflag_);
    PHISTTEST_MVEC_CREATE(&linearV_out_exact_,linearMap_,1,&iflag_);
    EXPECT_EQ(0,iflag_);

    PHISTTEST_MVEC_CREATE(&repartV_in_,repartMap_,1,&iflag_);
    EXPECT_EQ(0,iflag_);
    PHISTTEST_MVEC_CREATE(&repartV_out_,repartMap_,1,&iflag_);
    EXPECT_EQ(0,iflag_);
    PHISTTEST_MVEC_CREATE(&repartV_out_exact_,repartMap_,1,&iflag_);
    EXPECT_EQ(0,iflag_);

    // extract raw view for filling the input and precalculated result vector
    double *V_in_raw=NULL, *V_out_exact_raw=NULL;
    phist_lidx ldV_in, ldV_out;
    phist_Dmvec_extract_view(linearV_in_,&V_in_raw,&ldV_in,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_Dmvec_extract_view(linearV_out_exact_,&V_out_exact_raw,&ldV_out,&iflag_);
    EXPECT_EQ(0,iflag_);

#ifdef PHIST_MVECS_ROW_MAJOR
    EXPECT_EQ(1,ldV_in);
    EXPECT_EQ(1,ldV_out);
#endif

    // fill the input vector and compute exact result of A*x in the linear map
    phist_gidx ilower, iupper;
    phist_map_get_ilower(linearMap_,&ilower,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_map_get_iupper(linearMap_,&iupper,&iflag_);
    EXPECT_EQ(0,iflag_);

    phist_gidx DIM;
    ghost_lidx conf_spinZ[3] = {_L_,_L_/2,0};
    SpinChainSZ( -2, conf_spinZ, &DIM, NULL, NULL);

    matfuncs_info_t info;
    SpinChainSZ( -1, NULL, NULL, &info, NULL);

    double      vals[_N_];
    ghost_gidx  cols[_N_];
  
    for (ghost_gidx i=ilower; i<=iupper; i++)
    {
      ghost_lidx iloc = (ghost_lidx)(i-ilower);
      V_in_raw[iloc]=(double)(i+1);
      V_out_exact_raw[iloc]=0.0;
      
      ghost_lidx nnz;
      SpinChainSZ( i, &nnz, cols, (void*)vals, NULL);
      for (int j=0; j<nnz; j++)
      {
        V_out_exact_raw[iloc] += vals[j]*((double)(cols[j]+1)); 
      }
    }

    // clean up
    SpinChainSZ( -3, conf_spinZ, &DIM, NULL, NULL);
    
    // upload V_in and V_out_exact to device (if applicable)
    phist_Dmvec_to_device(linearV_in_,&iflag_);
    EXPECT_EQ(0,iflag_);

    phist_Dmvec_to_device(linearV_out_exact_,&iflag_);
    EXPECT_EQ(0,iflag_);
    
    // fill the remaining vectors with garbage, should all be overwritten by 
    // mvec_to_mvec and sparseMat_times_mvec in the tests below
    phist_Dmvec_put_value(linearV_out_,-23.0e23,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_Dmvec_put_value(repartV_in_,9999.0e99,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_Dmvec_put_value(repartV_out_,42.9e3,&iflag_);
    EXPECT_EQ(0,iflag_);
  }

  static void SetUpTestCase()
  {
    TestWithType<double>::SetUpTestCase();
    KernelTestWithMap<_N_>::SetUpTestCase();
  }

  /*! Clean up.
   */
  virtual void TearDown()
  {
    phist_DsparseMat_delete(linearA_,&iflag_);
    phist_Dmvec_delete(linearV_in_,&iflag_);
    phist_Dmvec_delete(linearV_out_,&iflag_);
    phist_Dmvec_delete(linearV_out_exact_,&iflag_);

    phist_DsparseMat_delete(repartA_,&iflag_);
    phist_Dmvec_delete(repartV_in_,&iflag_);
    phist_Dmvec_delete(repartV_out_,&iflag_);
    phist_Dmvec_delete(repartV_out_exact_,&iflag_);
  }

  static void TearDownTestCase()
  {
    KernelTestWithMap<_N_>::TearDownTestCase();
  }


protected:

  // two different maps: the regular one and the one obtained as row distribution
  // when asking for PHIST_SPARSEMAT_PERM_GLOBAL
  phist_const_map_ptr linearMap_ = NULL, repartMap_ = NULL;
  
  // input, result and exact result of spMVM in linear map
  phist_Dmvec_ptr linearV_in_ = NULL, linearV_out_ = NULL, linearV_out_exact_ = NULL;

  // input, result and exact result of spMVM in repartitioned map
  phist_Dmvec_ptr repartV_in_ = NULL, repartV_out_ = NULL, repartV_out_exact_ = NULL;
  
  // spinSZ[L] matrix with and without repartitioning
  phist_DsparseMat_ptr linearA_ = NULL, repartA_ = NULL;
  
};

  TEST_F(DSparseMatRepartTest,perm_symmetric)
  {
    // check that the domain and range map are the same for both the linear and repartitioned matrix
    // we do not impose anything on the row and column map since we don't care what happens internally
    // as long as the in- and output vectors have the same map.
    phist_const_map_ptr repart_range_map, linear_range_map;
    SUBR(sparseMat_get_range_map)(linearA_,&linear_range_map,&iflag_);
    ASSERT_EQ(0,iflag_);
    phist_maps_compatible(linearMap_,linear_range_map,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(sparseMat_get_range_map)(repartA_,&repart_range_map,&iflag_);
    phist_maps_compatible(repartMap_,repart_range_map,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(0,iflag_);
  }

  // tests that mvec_to_mvec recovers the exact vector from the permuted one.
  // if applied twice with exchanged arguments
  TEST_F(DSparseMatRepartTest, redist_vectors)
  {
    if (!typeImplemented_) return;
    phist_Dmvec_to_mvec(linearV_in_, repartV_in_,&iflag_);
    ASSERT_EQ(0,iflag_);
    // transfer  and unpermute the permuted vector to recover the original
    phist_Dmvec_to_mvec(repartV_in_, linearV_out_,&iflag_);    
#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"linearV_in_:\n");
    SUBR(mvec_print)(linearV_in_,&iflag_);
    ASSERT_EQ(0,iflag_);

    PHIST_SOUT(PHIST_DEBUG,"repartV_in_:\n");
    SUBR(mvec_print)(repartV_in_,&iflag_);
    ASSERT_EQ(0,iflag_);

    PHIST_SOUT(PHIST_DEBUG,"linearV_out_:\n");
    SUBR(mvec_print)(linearV_out_,&iflag_);
    ASSERT_EQ(0,iflag_);
#endif
    ASSERT_REAL_EQ(1.0,MvecsEqual(linearV_in_,linearV_out_));
  }

  TEST_F(DSparseMatRepartTest, standard_spmvm)
  {
    if (!typeImplemented_) return;
    phist_DsparseMat_times_mvec(1.0,linearA_,linearV_in_,0.0,linearV_out_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_REAL_EQ(1.0,MvecsEqual(linearV_out_,linearV_out_exact_));
  }


  TEST_F(DSparseMatRepartTest, permuted_spmvm)
  {
    if (!typeImplemented_) return;

    // first permute the input vector into the domain map of the matrix, which in this case
    // is the same as the row map of the repartitioned A.
    phist_Dmvec_to_mvec(linearV_in_,repartV_in_,&iflag_);
    ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"linearV_in_:\n");
    SUBR(mvec_print)(linearV_in_,&iflag_);
    ASSERT_EQ(0,iflag_);

    PHIST_SOUT(PHIST_DEBUG,"repartV_in_:\n");
    SUBR(mvec_print)(repartV_in_,&iflag_);
    ASSERT_EQ(0,iflag_);
#endif
    
    // perform spMVM with repartitioned/permuted matrix and vectors
    phist_DsparseMat_times_mvec(1.0,repartA_,repartV_in_,0.0,repartV_out_,&iflag_);
    ASSERT_EQ(0,iflag_);
    // unpermute the result
    phist_Dmvec_to_mvec(repartV_out_, linearV_out_, &iflag_);
    ASSERT_EQ(0,iflag_);

#if PHIST_OUTLEV>=PHIST_DEBUG
    PHIST_SOUT(PHIST_DEBUG,"repartV_out_:\n");
    SUBR(mvec_print)(repartV_out_,&iflag_);
    ASSERT_EQ(0,iflag_);
    PHIST_SOUT(PHIST_DEBUG,"linearV_out_:\n");
    SUBR(mvec_print)(linearV_out_,&iflag_);
    ASSERT_EQ(0,iflag_);
    PHIST_SOUT(PHIST_DEBUG,"linearV_out_exact_:\n");
    SUBR(mvec_print)(linearV_out_exact_,&iflag_);
    ASSERT_EQ(0,iflag_);
#endif

    // check correctness
    ASSERT_REAL_EQ(1.0,MvecsEqual(linearV_out_,linearV_out_exact_));    
  }

