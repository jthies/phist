#ifndef CLASSNAME
#error "file not included correctly."
#endif

// the case A=0 does not work
// for ghost because the binCRS format
// doesn't handle it correctly, it seems
#ifdef PHIST_KERNEL_LIB_GHOST
#ifndef SKIP_ZERO_MAT
#define SKIP_ZERO_MAT
#endif
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_>,
                 public KernelTestWithSdMats<_ST_,_NV_,_NV_>
{

  public:
  
  typedef KernelTestWithVectors<_ST_,_N_,_NV_> VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_> MTest;
  typedef KernelTestWithType< _MT_ > MT_Test;

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    VTest::SetUp();
    MTest::SetUp();
    
    if (typeImplemented_)
    {
      SUBR(read_mat)("speye",nglob_,&A1_,&iflag_);
#ifndef SKIP_ZERO_MAT
      SUBR(read_mat)("spzero",nglob_,&A0_,&iflag_);
#else
      A0_=A1_;
#endif
      SUBR(read_mat)("sprandn",nglob_,&A2_,&iflag_);
      SUBR(read_mat)("sprandn_nodiag",nglob_,&A3_,&iflag_);
#ifndef SKIP_ZERO_MAT
      SUBR(read_mat)("spshift",nglob_,&A4_,&iflag_);
#else
      A4_ = A1_;
#endif
      
      if (A0_==NULL || A1_==NULL || A2_==NULL || A3_==NULL || A4_==NULL)
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
    VTest::TearDown();
    MTest::TearDown();
    if (typeImplemented_)
    {
#ifndef SKIP_ZERO_MAT
      ASSERT_EQ(0,delete_mat(A0_));
#endif      
      ASSERT_EQ(0,delete_mat(A1_));
      ASSERT_EQ(0,delete_mat(A2_));
      ASSERT_EQ(0,delete_mat(A3_));
#ifndef SKIP_ZERO_MAT
      ASSERT_EQ(0,delete_mat(A4_));
#endif
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
    SUBR(mvec_create)(&vec1_,domain_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec1_,&vec1_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    SUBR(mvec_create)(&vec2_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec2_,&vec2_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    SUBR(mvec_create)(&vec3_,range_map,nvec_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_extract_view)(vec3_,&vec3_vp_,&lda,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(lda,lda_);

    phist_map_get_local_length(domain_map, &nloc_, &iflag_);
    ASSERT_EQ(0,iflag_);
  }
}


  void test_sparseMat_times_mvec_vadd_mvec(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ shifts[_NV_], _ST_ beta)
  {
    if( !typeImplemented_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts, vec1_, beta, vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    _ST_ alpha_shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      alpha_shifts[i] = alpha*shifts[i];

    SUBR(sparseMat_times_mvec)(alpha, A, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_vadd_mvec)(alpha_shifts, vec1_, st::one(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_,vflag_), sqrt(mt::eps()));
  }

#if _NV_ > 1
  void test_sparseMat_times_mvec_with_views(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ beta, int imin, int imax)
  {
    if( !typeImplemented_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    // create views
    TYPE(mvec_ptr) vec1_view = NULL;
    SUBR(mvec_view_block)(vec1_, &vec1_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);
    TYPE(mvec_ptr) vec2_view = NULL;
    SUBR(mvec_view_block)(vec2_, &vec2_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec)(alpha, A, vec1_view, beta, vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);

    // make sure nothing changed outside of viewed block
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,imin,lda_,stride_,vflag_));
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imax+1,lda_),vec3_vp_+VIDX(0,imax+1,lda_),nloc_,nvec_-imax-1,lda_,stride_,vflag_));

    // calculation for full block as reference
    SUBR(sparseMat_times_mvec)(alpha, A, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imin,lda_),vec3_vp_+VIDX(0,imin,lda_),nloc_,imax-imin+1,lda_,stride_,vflag_), sqrt(mt::eps()));

    // delete view
    SUBR(mvec_delete)(vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_delete)(vec1_view, &iflag_);
    ASSERT_EQ(0, iflag_);
  }

  void test_sparseMat_times_mvec_vadd_mvec_with_views(_ST_ alpha, TYPE(const_sparseMat_ptr) A, _ST_ shifts[_NV_], _ST_ beta, int imin, int imax)
  {
    if( !typeImplemented_ )
      return;

    // set up mvecs
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(vec2_, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_add_mvec)(st::one(), vec2_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    // create a view
    TYPE(mvec_ptr) vec1_view = NULL;
    SUBR(mvec_view_block)(vec1_, &vec1_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);
    TYPE(mvec_ptr) vec2_view = NULL;
    SUBR(mvec_view_block)(vec2_, &vec2_view, imin, imax, &iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts+imin, vec1_view, beta, vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);

    // make sure nothing changed outside of viewed block
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_,vec3_vp_,nloc_,imin,lda_,stride_,vflag_));
    ASSERT_REAL_EQ(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imax+1,lda_),vec3_vp_+VIDX(0,imax+1,lda_),nloc_,nvec_-imax-1,lda_,stride_,vflag_));

    // calculation for full block as reference
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A, shifts, vec1_, beta, vec3_, &iflag_);
    ASSERT_EQ(0, iflag_);

    ASSERT_NEAR(mt::one(), ArraysEqual(vec2_vp_+VIDX(0,imin,lda_),vec3_vp_+VIDX(0,imin,lda_),nloc_,imax-imin+1,lda_,stride_,vflag_), sqrt(mt::eps()));

    // delete view
    SUBR(mvec_delete)(vec2_view, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_delete)(vec1_view, &iflag_);
    ASSERT_EQ(0, iflag_);
  }
#endif

  void test_carp_kernel(TYPE(const_sparseMat_ptr) A)
  {
    if (!typeImplemented_) return;
    // ColPack has trouble with the tiny
    // local matrices that occur when partitioning a
    // 25x25 matrix (TODO: larger CrsMat tests)
    if (mpi_size_>1) return;

    // setup the CARP kernel and get the required data structures:
    MT* nrms_ai2i=NULL;
    void* aux=NULL;
    int numSys=1;
    MT sigma_r=0.0;
    MT sigma_i=0.0;
    MT omega=1.0;
    SUBR(carp_setup)(A,numSys,&sigma_r,&sigma_i,
        &nrms_ai2i, &aux, &iflag_);
    if (iflag_==-99) return; // CARP not implemented
    ASSERT_EQ(0,iflag_);
    
    // perform single CARP forward/backward sweep with B=0 and X=v2.
    
    // set up mvecs
    TYPE(mvec_ptr) B=vec1_;
    TYPE(mvec_ptr) Xr=vec2_;
    TYPE(mvec_ptr) Xi=vec3_;
    SUBR(mvec_put_value)(B, st::zero(),&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(Xr, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_random)(Xi, &iflag_);
    ASSERT_EQ(0, iflag_);
    
    MT norms_X0[_NV_];
    MT norms_X1[_NV_];
    MT norms_X2[_NV_];

    // backup X
    TYPE(mvec_ptr) Xr_bak=NULL;
    TYPE(mvec_ptr) Xi_bak=NULL;
    SUBR(mvec_create)(&Xr_bak,map_,_NV_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_create)(&Xi_bak,map_,_NV_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(), Xr, st::zero(), Xr_bak, &iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_add_mvec)(st::one(), Xi, st::zero(), Xi_bak, &iflag_);
    ASSERT_EQ(0, iflag_);
    
    SUBR(carp_sweep)(A, 1, &sigma_r, &sigma_i,B,&Xr,&Xi,
          nrms_ai2i,aux,&omega,&iflag_);
    ASSERT_EQ(0, iflag_);
    
    // compute norms before and after the sweep
    SUBR(mvec_norm2)(Xr_bak,norms_X0,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(mvec_norm2)(Xr,norms_X1,&iflag_);
    ASSERT_EQ(0, iflag_);

    PHIST_SOUT(PHIST_VERBOSE,"||X||_2 before and after CARP sweep:\n");
    for (int i=0; i<nvec_; i++) PHIST_SOUT(PHIST_VERBOSE,
        "%12.8e\t%12.8e\n", norms_X0[i], norms_X1[i]);

    // TEST 1: row projection must be non-increasing on the vectors norm.
    for (int i=0; i<_NV_; i++)
    {
      ASSERT_TRUE(norms_X0[i]>=norms_X1[i]);
    }
    

    // Check that it works if B=NULL is passed in
    SUBR(mvec_add_mvec)(st::one(),Xr_bak,st::zero(),Xr,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(carp_sweep)(A, 1, &sigma_r, &sigma_i,NULL,&Xr,&Xi,
          nrms_ai2i,aux,&omega,&iflag_);
    ASSERT_EQ(0, iflag_);

    SUBR(mvec_norm2)(Xr,norms_X2,&iflag_);
    ASSERT_EQ(0, iflag_);

    PHIST_SOUT(PHIST_VERBOSE,"||X||_2 after CARP with B=0 and B=NULL:\n");
    for (int i=0; i<nvec_; i++) PHIST_SOUT(PHIST_VERBOSE,
        "%12.8e\t%12.8e\n", norms_X1[i], norms_X2[i]);
    
    ASSERT_REAL_EQ(mt::one(),MT_Test::ArraysEqual(norms_X1,norms_X2,nvec_,1,nvec_,1,vflag_));
      
    // test that the forward/backward CARP operator is symmetric
    SUBR(mvecT_times_mvec)(st::one(),Xr,Xr_bak,st::zero(),mat1_,&iflag_);
    ASSERT_EQ(0, iflag_);
    SUBR(sdMat_from_device)(mat1_,&iflag_);
    ASSERT_EQ(0, iflag_);      
    MT nonsymm=mt::zero();
    for (int i=1; i<nvec_; i++)
    {
      for (int j=1; j<nvec_; j++)
      {
        nonsymm = std::max(nonsymm,
                st::abs(mat1_vp_[i*lda_+j]-mat1_vp_[j*lda_+i]));
      }
    }
    ASSERT_NEAR(mt::one(),mt::one()+nonsymm,100*mt::eps());

    SUBR(carp_destroy)(A,1,nrms_ai2i,aux,&iflag_);
    ASSERT_EQ(0, iflag_);

  SUBR(mvec_delete)(Xr_bak,&iflag_);
    ASSERT_EQ(0,iflag_);

  SUBR(mvec_delete)(Xi_bak,&iflag_);
    ASSERT_EQ(0,iflag_);

  }

TYPE(sparseMat_ptr) A0_; // all zero matrix
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

_MT_ const_row_sum_test(TYPE(sparseMat_ptr) A)
  {
    if (typeImplemented_ && haveMats_)
    {
      _ST_ val = random_number();
      global_sum(&val,1,mpi_comm_);
      SUBR(mvec_put_value)(vec1_,val,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A,vec1_,st::zero(),vec2_,&iflag_);
      if (iflag_) return (_MT_)iflag_;
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec1_,&iflag_);
      SUBR(mvec_print)(vec2_,&iflag_);
#endif
      return ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,val,vflag_);
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

#ifndef SKIP_ZERO_MAT
  TEST_F(CLASSNAME, A0_times_mvec) 
    {
    if (typeImplemented_ && haveMats_)
      {
      // matrices may have different maps
      rebuildVectors(A0_);

      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A0_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,0.0,vflag_));
      }
    }
#endif

  TEST_F(CLASSNAME, A1_times_mvec)
    {
    if (typeImplemented_ && haveMats_)
      {
      // matrices may have different maps
      rebuildVectors(A1_);

      ST alpha, beta;
      //I*X=X?
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(st::one(),A1_,vec1_,st::zero(),vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_,vflag_));

      //alpha*I*X=alpha*X?
      alpha = random_number();
      beta=st::zero();
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      SUBR(sparseMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_scale)(vec1_,alpha,&iflag_);
      ASSERT_EQ(0,iflag_);
      ASSERT_NEAR(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_,vflag_),100*mt::one());

      //0*I*X+beta*Y = beta*Y? 
      alpha=st::zero(); 
      beta=random_number();
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
#endif
      SUBR(mvec_random)(vec1_,&iflag_); 
      SUBR(mvec_random)(vec2_,&iflag_); 
#if PHIST_OUTLEV>=PHIST_DEBUG
      std::cout << "input="<<std::endl;
      SUBR(mvec_print)(vec1_,&iflag_);
      std::cout << "output, before="<<std::endl;
      SUBR(mvec_print)(vec2_,&iflag_);
#endif
      // v3=beta*v2 
      SUBR(mvec_add_mvec)(beta,vec2_,st::zero(),vec3_,&iflag_); 
      ASSERT_EQ(0,iflag_); 
      // v2 = 0*v1 + beta*v2 (=v3) 
      SUBR(sparseMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&iflag_); 
      ASSERT_EQ(0,iflag_); 
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&iflag_);
      SUBR(mvec_print)(vec3_,&iflag_);
#endif
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_,vflag_));

      //I*X+beta*Y = X+beta*Y?
      alpha = st::one();
      beta = random_number();
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
#endif
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      std::cout << "input="<<std::endl;
      SUBR(mvec_print)(vec1_,&iflag_);
      std::cout << "output, before="<<std::endl;
      SUBR(mvec_print)(vec2_,&iflag_);
#endif
      //v3=v1+beta*v2
      SUBR(mvec_to_mvec)(vec1_,vec3_,&iflag_);
      if( iflag_ == PHIST_NOT_IMPLEMENTED )
      {
        SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec3_,&iflag_);
      }
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_scale)(vec3_, alpha, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // v2 = v1 + beta*v2 (=alpha*v1+v3)
      SUBR(sparseMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&iflag_);
      SUBR(mvec_print)(vec3_,&iflag_);
#endif
      ASSERT_NEAR(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_,vflag_),100*mt::eps());

      //alpha*I*X+beta*Y = alpha*X+beta*Y?
      alpha = random_number();
      beta = random_number();
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A=I, alpha="<<alpha<<", beta="<<beta<<std::endl;
#endif
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      std::cout << "input="<<std::endl;
      SUBR(mvec_print)(vec1_,&iflag_);
      std::cout << "output, before="<<std::endl;
      SUBR(mvec_print)(vec2_,&iflag_);
#endif
       // v3=alpha*v1+beta*v2
      SUBR(mvec_to_mvec)(vec1_,vec3_,&iflag_);
      if( iflag_ == PHIST_NOT_IMPLEMENTED )
      {
        SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec3_,&iflag_);
      }
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_scale)(vec3_, alpha, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(beta,vec2_,st::one(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // v2 = alpha*v1 + beta*v2 (=alpha*v1+v3)
      SUBR(sparseMat_times_mvec)(alpha,A1_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&iflag_);
      SUBR(mvec_print)(vec3_,&iflag_);
#endif
      ASSERT_NEAR(mt::one(),ArraysEqual(vec2_vp_,vec3_vp_,nloc_,nvec_,lda_,stride_,vflag_),100*mt::eps());
      }
    }

#if(_NV_>1)
  TEST_F(CLASSNAME, A1_times_mvec_using_two_views_of_the_same_vec)
  {
    if (typeImplemented_ && A1_!=NULL)
    {
      // matrices may have different maps
      rebuildVectors(A1_);

      ST alpha, beta;
      //I*X=X?
      SUBR(mvec_random)(vec1_,&iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      
      TYPE(mvec_ptr) v_in=NULL, v_out=NULL;
      lidx_t nv =  _NV_/2;
      lidx_t offs = _NV_%2;
      SUBR(mvec_view_block)(vec1_,&v_in,offs,offs+nv-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_view_block)(vec1_,&v_out,offs+nv,offs+2*nv-1,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      _ST_ *v_in_vp, *v_out_vp;

// we could extract the pointers from the view, but this way we
// have more control and make sure the data is where it should be:
     v_in_vp = vec1_vp_+VIDX(0,offs,lda_);
     v_out_vp = vec1_vp_+VIDX(0,offs+nv,lda_);
#if(PHIST_OUTLEV>=PHIST_DEBUG)
      PHIST_DEB("all random:\n");
      SUBR(mvec_print)(vec1_,&iflag_);
#endif
      SUBR(sparseMat_times_mvec)(st::one(),A1_,v_in,st::zero(),v_out,&iflag_);
      ASSERT_EQ(0,iflag_);
#if(PHIST_OUTLEV>=PHIST_DEBUG)
      PHIST_DEB("with two identical blocks [%d..%d] and [%d..%d]:\n",
        offs, offs+nv-1,offs+nv,offs+2*nv-1);
      SUBR(mvec_print)(vec1_,&iflag_);
#endif
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(v_in_vp,v_out_vp,nloc_,nv,lda_,stride_,vflag_));
    }
    else
    {
      PHIST_SOUT(PHIST_INFO,"test skipped\n");
    }
  }
#endif

  TEST_F(CLASSNAME, A2_times_mvec)
    {
      // matrices may have different maps
      rebuildVectors(A2_);

    // we allow a tolerance here because the matrices may have errors in the
    // last digit and we can't get the test to pass otherwise.
    ASSERT_NEAR(mt::one(),const_row_sum_test(A2_),100*mt::eps());
    }

  TEST_F(CLASSNAME, A3_times_mvec)
    {
    // matrices may have different maps
    rebuildVectors(A3_);

    // we allow a tolerance here because the matrices may have errors in the
    // last digit and we can't get the test to pass otherwise.
    ASSERT_NEAR(mt::one(),const_row_sum_test(A3_),100*mt::eps());
    }

#ifndef SKIP_ZERO_MAT
  TEST_F(CLASSNAME, shift_mvec)
  {
    if( typeImplemented_ )
    {
      // matrices may have different maps
      rebuildVectors(A4_);

      _ST_ alpha = st::one();
      _ST_ beta = st::zero();

      // create new map (which has correct order!)
      map_ptr_t map;
      phist_map_create(&map, comm_, nglob_, &iflag_);
      ASSERT_EQ(0,iflag_);
      lidx_t nloc = 0;
      phist_map_get_local_length(map,&nloc,&iflag_);
      ASSERT_EQ(0,iflag_);
      gidx_t ilower = 0;
      phist_map_get_ilower(map,&ilower,&iflag_);
      ASSERT_EQ(0,iflag_);

      // create a vector with this map
      TYPE(mvec_ptr) orderedVec = NULL;
      SUBR(mvec_create)(&orderedVec, map, nvec_, &iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ *orderedVec_vp = NULL;
      lidx_t lda = 0;
      SUBR(mvec_extract_view)(orderedVec, &orderedVec_vp, &lda, &iflag_);
      ASSERT_EQ(0,iflag_);

      // setup recognizable input
      for(int i = 0; i < nloc; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          orderedVec_vp[VIDX(i,j,lda)] = (_ST_)(ilower+i + j*nglob_);
        }
      }
      // copy to vec1_
      SUBR(mvec_to_mvec)(orderedVec, vec1_, &iflag_);
      if( iflag_ == PHIST_NOT_IMPLEMENTED )
      {
        SUBR(mvec_add_mvec)(st::one(),orderedVec,st::zero(),vec1_,&iflag_);
      }
      ASSERT_EQ(0,iflag_);


      // apply our shift matrix
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A='shift', alpha=1, beta=0"<<std::endl;
#endif
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      SUBR(sparseMat_times_mvec)(alpha,A4_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      // copy to orderedVec
      SUBR(mvec_to_mvec)(vec2_, orderedVec, &iflag_);
      if( iflag_ == PHIST_NOT_IMPLEMENTED )
      {
        SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),orderedVec,&iflag_);
      }
      ASSERT_EQ(0,iflag_);


      // check result
      for(int i = 0; i < nglob_; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          if( i >= ilower && i < ilower+nloc )
          {
            ASSERT_REAL_EQ((i+1)%nglob_ + j*nglob_,st::real(orderedVec_vp[VIDX(i-ilower,j,lda)]));
          }
          else
          {
            // also assert elements of other processes, s.t. all processes do the same number of asserts!
            ASSERT_REAL_EQ((i+1)%nglob_ + j*nglob_, (i+1)%nglob_ + j*nglob_);
          }
        }
      }
    }
  }
#endif


#if defined(PHIST_KERNEL_LIB_GHOST)&&defined(PHIST_USE_SELL)&&(PHIST_SELL_SIGMA>1)
  // this test is broken with local reordering because we can't create
  // an "unpermuted" vector in the map of A with the present interface.
  TEST_F(CLASSNAME, DISABLED_A2_precalc_result)
#else
  TEST_F(CLASSNAME, A2_precalc_result)
#endif
  {
    if( typeImplemented_ )
    {
      // matrices may have different maps
      rebuildVectors(A2_);

      _ST_ alpha = st::one();
      _ST_ beta = st::zero();

      // create new map (which has correct order!)
      map_ptr_t map;
      phist_map_create(&map, comm_, nglob_, &iflag_);
      ASSERT_EQ(0,iflag_);
      lidx_t nloc = 0;
      phist_map_get_local_length(map,&nloc,&iflag_);
      ASSERT_EQ(0,iflag_);
      gidx_t ilower = 0;
      phist_map_get_ilower(map,&ilower,&iflag_);
      ASSERT_EQ(0,iflag_);

      // create a vector with this map
      TYPE(mvec_ptr) orderedVec = NULL;
      SUBR(mvec_create)(&orderedVec, map, nvec_, &iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ *orderedVec_vp = NULL;
      lidx_t lda = 0;
      SUBR(mvec_extract_view)(orderedVec, &orderedVec_vp, &lda, &iflag_);
      ASSERT_EQ(0,iflag_);

      // setup recognizable input
      for(int i = 0; i < nloc; i++)
      {
        for(int j = 0; j < nvec_; j++)
        {
          orderedVec_vp[VIDX(i,j,lda)] = (_ST_)(ilower+i + j*nglob_);
        }
      }
      // copy to vec1_
      SUBR(mvec_to_mvec)(orderedVec, vec1_, &iflag_);
      if( iflag_ == PHIST_NOT_IMPLEMENTED )
      {
        SUBR(mvec_add_mvec)(st::one(),orderedVec,st::zero(),vec1_,&iflag_);
      }
      ASSERT_EQ(0,iflag_);


      // apply our shift matrix
#if PHIST_OUTLEV>=PHIST_INFO
      std::cout << "MVM with A='rand', alpha=1, beta=0"<<std::endl;
#endif
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif
      SUBR(sparseMat_times_mvec)(alpha,A2_,vec1_,beta,vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      SUBR(mvec_print)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
#endif
 
#if _N_ == 25 && _NV_ == 1
#ifdef IS_COMPLEX
      _ST_ precalc_result[_N_*_NV_] = {
        _ST_(1.072230172746e+01,-1.031665649804e+00), 
        _ST_(2.131001265990e+01,-3.282812619721e+00), 
        _ST_(2.015641860460e+01,-1.541996705535e+01), 
        _ST_(2.171736316821e+01,6.518151979141e+00), 
        _ST_(2.615219127525e+00,-3.620215930236e+00), 
        _ST_(1.228901252787e+01,-4.528370474365e+00), 
        _ST_(-6.543281193088e+00,4.169322949310e+00), 
        _ST_(1.513919626135e+01,-3.968740334506e+00), 
        _ST_(3.010875585729e+01,-7.275935201946e+00), 
        _ST_(1.856357361375e+01,-6.082795919035e-01), 
        _ST_(4.884622853589e+00,-1.246056153190e+01), 
        _ST_(9.724002250734e+00,1.854295108448e+00), 
        _ST_(1.295234521942e+01,4.797517156124e+00), 
        _ST_(1.922176454031e+01,-6.554985875651e+00), 
        _ST_(9.937028268565e+00,-3.048436869856e+00), 
        _ST_(8.417823023929e+00,1.504269262034e+00), 
        _ST_(2.827516113199e+00,-6.272713645344e+00), 
        _ST_(1.542580758006e+01,-5.650405298966e+00), 
        _ST_(1.507749351511e+01,-4.585936779631e+00), 
        _ST_(1.644503601436e+01,-6.053982924770e+00), 
        _ST_(3.407669499649e+01,1.638938366859e+01), 
        _ST_(2.078494495946e+01,5.856000872126e-01), 
        _ST_(1.420728160837e+01,1.050663583365e+00), 
        _ST_(1.251387547773e+01,-4.875170138989e-01), 
        _ST_(3.325495196562e+00,-1.806246226644e+01) };
#else
      _ST_ precalc_result[_N_*_NV_] = {
        -0.12861779756316816,
        16.076787993876188,
        20.072300560349166,
        17.898973758721571,
        19.186951140368311,
        -52.388079187547561,
        12.117945192440743,
        11.417408852483819,
        18.080795442101312,
        13.972038117529888,
        4.3052769861597753,
        12.535799082066269,
        5.9149777593318920,
        8.3650228590087750E-002,
        10.824168487154179,
        44.098112302021931,
        80.769125851754040,
        10.243308300563234,
        -40.240733293758780,
        16.598578972374874,
        18.400975207367228,
        36.482850494569007,
        13.997556100937755,
        -8.1932230989204431,
        9.6675923086516082 };
#endif
#elif _N_ == 25 && _NV_ == 4
#ifdef IS_COMPLEX
      _ST_ precalc_result[_N_*_NV_] = {
        _ST_(1.072230172746e+01,-1.031665649804e+00), _ST_(3.572230172746e+01,-1.031665649804e+00), _ST_(6.072230172746e+01,-1.031665649804e+00), _ST_(8.572230172746e+01,-1.031665649804e+00),
        _ST_(2.131001265990e+01,-3.282812619721e+00), _ST_(4.631001265990e+01,-3.282812619721e+00), _ST_(7.131001265990e+01,-3.282812619721e+00), _ST_(9.631001265990e+01,-3.282812619721e+00),
        _ST_(2.015641860460e+01,-1.541996705535e+01), _ST_(4.515641860460e+01,-1.541996705535e+01), _ST_(7.015641860460e+01,-1.541996705535e+01), _ST_(9.515641860460e+01,-1.541996705535e+01),
        _ST_(2.171736316821e+01,6.518151979141e+00), _ST_(4.671736316821e+01,6.518151979141e+00), _ST_(7.171736316821e+01,6.518151979141e+00), _ST_(9.671736316821e+01,6.518151979141e+00),
        _ST_(2.615219127525e+00,-3.620215930236e+00), _ST_(2.761521912752e+01,-3.620215930236e+00), _ST_(5.261521912752e+01,-3.620215930236e+00), _ST_(7.761521912752e+01,-3.620215930236e+00),
        _ST_(1.228901252787e+01,-4.528370474365e+00), _ST_(3.728901252787e+01,-4.528370474365e+00), _ST_(6.228901252787e+01,-4.528370474365e+00), _ST_(8.728901252787e+01,-4.528370474365e+00),
        _ST_(-6.543281193088e+00,4.169322949310e+00), _ST_(1.845671880691e+01,4.169322949310e+00), _ST_(4.345671880691e+01,4.169322949310e+00), _ST_(6.845671880691e+01,4.169322949310e+00),
        _ST_(1.513919626135e+01,-3.968740334506e+00), _ST_(4.013919626135e+01,-3.968740334506e+00), _ST_(6.513919626135e+01,-3.968740334506e+00), _ST_(9.013919626135e+01,-3.968740334506e+00),
        _ST_(3.010875585729e+01,-7.275935201946e+00), _ST_(5.510875585729e+01,-7.275935201946e+00), _ST_(8.010875585729e+01,-7.275935201946e+00), _ST_(1.051087558573e+02,-7.275935201946e+00),
        _ST_(1.856357361375e+01,-6.082795919035e-01), _ST_(4.356357361375e+01,-6.082795919035e-01), _ST_(6.856357361375e+01,-6.082795919035e-01), _ST_(9.356357361375e+01,-6.082795919035e-01),
        _ST_(4.884622853589e+00,-1.246056153190e+01), _ST_(2.988462285359e+01,-1.246056153190e+01), _ST_(5.488462285359e+01,-1.246056153190e+01), _ST_(7.988462285359e+01,-1.246056153190e+01),
        _ST_(9.724002250734e+00,1.854295108448e+00), _ST_(3.472400225073e+01,1.854295108448e+00), _ST_(5.972400225073e+01,1.854295108448e+00), _ST_(8.472400225073e+01,1.854295108448e+00),
        _ST_(1.295234521942e+01,4.797517156124e+00), _ST_(3.795234521942e+01,4.797517156124e+00), _ST_(6.295234521942e+01,4.797517156124e+00), _ST_(8.795234521942e+01,4.797517156124e+00),
        _ST_(1.922176454031e+01,-6.554985875651e+00), _ST_(4.422176454031e+01,-6.554985875651e+00), _ST_(6.922176454031e+01,-6.554985875651e+00), _ST_(9.422176454031e+01,-6.554985875651e+00),
        _ST_(9.937028268565e+00,-3.048436869856e+00), _ST_(3.493702826857e+01,-3.048436869856e+00), _ST_(5.993702826857e+01,-3.048436869856e+00), _ST_(8.493702826857e+01,-3.048436869856e+00),
        _ST_(8.417823023929e+00,1.504269262034e+00), _ST_(3.341782302393e+01,1.504269262034e+00), _ST_(5.841782302393e+01,1.504269262034e+00), _ST_(8.341782302393e+01,1.504269262034e+00),
        _ST_(2.827516113199e+00,-6.272713645344e+00), _ST_(2.782751611320e+01,-6.272713645344e+00), _ST_(5.282751611320e+01,-6.272713645344e+00), _ST_(7.782751611320e+01,-6.272713645344e+00),
        _ST_(1.542580758006e+01,-5.650405298966e+00), _ST_(4.042580758006e+01,-5.650405298966e+00), _ST_(6.542580758006e+01,-5.650405298966e+00), _ST_(9.042580758006e+01,-5.650405298966e+00),
        _ST_(1.507749351511e+01,-4.585936779631e+00), _ST_(4.007749351511e+01,-4.585936779631e+00), _ST_(6.507749351511e+01,-4.585936779631e+00), _ST_(9.007749351511e+01,-4.585936779631e+00),
        _ST_(1.644503601436e+01,-6.053982924770e+00), _ST_(4.144503601436e+01,-6.053982924770e+00), _ST_(6.644503601436e+01,-6.053982924770e+00), _ST_(9.144503601436e+01,-6.053982924770e+00),
        _ST_(3.407669499649e+01,1.638938366859e+01), _ST_(5.907669499649e+01,1.638938366859e+01), _ST_(8.407669499649e+01,1.638938366859e+01), _ST_(1.090766949965e+02,1.638938366859e+01),
        _ST_(2.078494495946e+01,5.856000872126e-01), _ST_(4.578494495946e+01,5.856000872126e-01), _ST_(7.078494495946e+01,5.856000872126e-01), _ST_(9.578494495946e+01,5.856000872125e-01),
        _ST_(1.420728160837e+01,1.050663583365e+00), _ST_(3.920728160837e+01,1.050663583365e+00), _ST_(6.420728160837e+01,1.050663583365e+00), _ST_(8.920728160837e+01,1.050663583365e+00),
        _ST_(1.251387547773e+01,-4.875170138989e-01), _ST_(3.751387547773e+01,-4.875170138989e-01), _ST_(6.251387547773e+01,-4.875170138989e-01), _ST_(8.751387547773e+01,-4.875170138989e-01),
        _ST_(3.325495196562e+00,-1.806246226644e+01), _ST_(2.832549519656e+01,-1.806246226644e+01), _ST_(5.332549519656e+01,-1.806246226644e+01), _ST_(7.832549519656e+01,-1.806246226644e+01) };
#else
      _ST_ precalc_result[_N_*_NV_] = {
        -0.12861779756316816, 24.871382202436834, 49.871382202436841, 74.871382202436862,
        16.076787993876188, 41.076787993876181, 66.076787993876167, 91.076787993876181,
        20.072300560349166, 45.072300560349177, 70.072300560349163, 95.072300560349149,
        17.898973758721571, 42.898973758721567, 67.898973758721596, 92.898973758721581,
        19.186951140368311, 44.186951140368315, 69.186951140368308, 94.186951140368308,
        -52.388079187547561, -27.388079187547561, -2.3880791875475609, 22.611920812452411,
        12.117945192440743, 37.117945192440750, 62.117945192440757, 87.117945192440772,
        11.417408852483819, 36.417408852483824, 61.417408852483824, 86.417408852483831,
        18.080795442101312, 43.080795442101319, 68.080795442101319, 93.080795442101319,
        13.972038117529888, 38.972038117529891, 63.972038117529905, 88.972038117529905,
        4.3052769861597753, 29.305276986159761, 54.305276986159726, 79.305276986159711,
        12.535799082066269, 37.535799082066255, 62.535799082066234, 87.535799082066234,
        5.9149777593318920, 30.914977759331908, 55.914977759331904, 80.914977759331919,
        8.3650228590087750E-002, 25.083650228590084, 50.083650228590088, 75.083650228590102,
        10.824168487154179, 35.824168487154182, 60.824168487154182, 85.824168487154196,
        44.098112302021931, 69.098112302021960, 94.098112302021960, 119.09811230202197,
        80.769125851754040, 105.76912585175401, 130.76912585175404, 155.76912585175393,
        10.243308300563234, 35.243308300563228, 60.243308300563228, 85.243308300563228,
        -40.240733293758780, -15.240733293758730, 9.7592667062412914, 34.759266706241320,
        16.598578972374874, 41.598578972374867, 66.598578972374881, 91.598578972374881,
        18.400975207367228, 43.400975207367232, 68.400975207367239, 93.400975207367239,
        36.482850494569007, 61.482850494569007, 86.482850494569036, 111.48285049456905,
        13.997556100937755, 38.997556100937771, 63.997556100937771, 88.997556100937743,
        -8.1932230989204431, 16.806776901079530, 41.806776901079516, 66.806776901079530,
        9.6675923086516082, 34.667592308651606, 59.667592308651606, 84.667592308651621 };
#endif
#else
  PHIST_SOUT(PHIST_INFO,"skipping test, no data available\n");
  ST* precalc_result=NULL;
  return;
#endif
      // copy to orderedVec
      SUBR(mvec_to_mvec)(vec2_, orderedVec, &iflag_);
      if( iflag_ == PHIST_NOT_IMPLEMENTED )
      {
        SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),orderedVec,&iflag_);
      }
      ASSERT_EQ(0,iflag_);


      // check result
      _MT_ err = 0;
      for(int i = 0; i < nloc; i++)
        for(int j = 0; j < nvec_; j++)
          err = std::max(err, st::abs( precalc_result[(ilower+i)*nvec_+j] - orderedVec_vp[VIDX(i,j,lda)] ));
      ASSERT_NEAR(err, mt::zero(), mt::sqrt(mt::eps()));
    }
  }


  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_only_scale)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::zero();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A2_, shifts, beta);
  }

#ifndef SKIP_ZERO_MAT
  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_only_vadd)
  {
    // matrices may have different maps
    rebuildVectors(A0_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::one();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A0_, shifts, beta);
  }
#endif

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_only_spmvm)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::zero();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A2_, shifts, beta);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec(alpha, A2_, shifts, beta);
  }

#if _NV_ > 1
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_0_0)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A2_, beta, 0, 0);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_0_0)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A2_, shifts, beta, 0, 0);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_1_1)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A2_, beta, 1, 1);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_1_1)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A2_, shifts, beta, 1, 1);
  }
#endif


#if _NV_ > 3
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_1_2)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A2_, beta, 1, 2);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_1_2)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A2_, shifts, beta, 1, 2);
  }
#endif


#if _NV_ > 3
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_2_3)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A2_, beta, 2, 3);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_2_3)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A2_, shifts, beta, 2, 3);
  }
#endif


#if _NV_ > 4
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_view_1_4)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_with_views(alpha, A2_, beta, 1, 4);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_view_1_4)
  {
    // matrices may have different maps
    rebuildVectors(A2_);

    _ST_ shifts[_NV_];
    for(int i = 0; i < _NV_; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();

    test_sparseMat_times_mvec_vadd_mvec_with_views(alpha, A2_, shifts, beta, 1, 4);
  }
#endif


#if _NV_ >= 4
  TEST_F(CLASSNAME, sparseMat_times_mvec_random_with_same_vec_views_0_1__2_3)
  {
    if( !typeImplemented_ )
      return;

    // matrices may have different maps
    rebuildVectors(A2_);

    // random data
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0,iflag_);

    // safe vec1_
    SUBR(mvec_add_mvec)(st::one(), vec1_, st::zero(), vec2_, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_to_mvec)(vec1_,vec3_,&iflag_);
    if( iflag_ == PHIST_NOT_IMPLEMENTED )
    {
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec3_,&iflag_);
    }
    ASSERT_EQ(0,iflag_);
    ASSERT_EQ(0,iflag_);

    // create views
    TYPE(mvec_ptr) vin = NULL;
    SUBR(mvec_view_block)(vec1_, &vin, 2, 3, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vout = NULL;
    SUBR(mvec_view_block)(vec1_, &vout, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vref = NULL;
    SUBR(mvec_view_block)(vec2_, &vref, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);

    // first generate reference data (safe calculation)
    SUBR(sparseMat_times_mvec)(alpha, A2_, vin, beta, vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    // calculation (unsafe aliasing!)
    SUBR(sparseMat_times_mvec)(alpha, A2_, vin, beta, vout, &iflag_);
    ASSERT_EQ(0,iflag_);

    // check vin
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_+VIDX(0,2,lda_),vec3_vp_+VIDX(0,2,lda_),nloc_,2,lda_,stride_,vflag_), mt::eps());
    // check result
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_,vec2_vp_,nloc_,2,lda_,stride_,vflag_), sqrt(mt::eps()));

    // delete views
    SUBR(mvec_delete)(vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vout, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vin, &iflag_);
    ASSERT_EQ(0,iflag_);
  }

  TEST_F(CLASSNAME, sparseMat_times_mvec_vadd_mvec_random_with_same_vec_views_0_1__2_3)
  {
    if( !typeImplemented_ )
      return;

    // matrices may have different maps
    rebuildVectors(A2_);

    // random data
    _ST_ shifts[2];
    for(int i = 0; i < 2; i++)
      shifts[i] = st::prand();
    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
    SUBR(mvec_random)(vec1_, &iflag_);
    ASSERT_EQ(0,iflag_);

    // safe vec1_
    SUBR(mvec_add_mvec)(st::one(), vec1_, st::zero(), vec2_, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(), vec1_, st::zero(), vec3_, &iflag_);
    ASSERT_EQ(0,iflag_);
    // create views
    TYPE(mvec_ptr) vin = NULL;
    SUBR(mvec_view_block)(vec1_, &vin, 2, 3, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vout = NULL;
    SUBR(mvec_view_block)(vec1_, &vout, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);
    TYPE(mvec_ptr) vref = NULL;
    SUBR(mvec_view_block)(vec2_, &vref, 0, 1, &iflag_);
    ASSERT_EQ(0,iflag_);

    // first generate reference data (safe calculation)
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A2_, shifts, vin, beta, vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    // calculation (unsafe aliasing!)
    SUBR(sparseMat_times_mvec_vadd_mvec)(alpha, A2_, shifts, vin, beta, vout, &iflag_);
    ASSERT_EQ(0,iflag_);

    // check vin
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_+VIDX(0,2,lda_),vec3_vp_+VIDX(0,2,lda_),nloc_,2,lda_,stride_,vflag_), mt::eps());
    // check result
    ASSERT_NEAR(mt::one(), ArraysEqual(vec1_vp_,vec2_vp_,nloc_,2,lda_,stride_,vflag_), sqrt(mt::eps()));

    // delete views
    SUBR(mvec_delete)(vref, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vout, &iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_delete)(vin, &iflag_);
    ASSERT_EQ(0,iflag_);
  }
#endif

//////////////////////////////////////////////////
// basic tests for the CARP kernel.             //
//////////////////////////////////////////////////
// tests disabled for now, they run on Emmy but not in the Jenkins builds at DLR
TEST_F(CLASSNAME, DISABLED_Carp_A1)
{
  test_carp_kernel(A1_);
}

TEST_F(CLASSNAME, DISABLED_Carp_A2)
{
  test_carp_kernel(A2_);
}

TEST_F(CLASSNAME, DISABLED_Carp_A3)
{
  test_carp_kernel(A3_);
}

