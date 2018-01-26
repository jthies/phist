/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,5>, 
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_,0>
  {

public:

  typedef KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,5> VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_,0> MTest;

#ifdef PHIST_HAVE_BELOS
  typedef Belos::MultiVecTraits< _ST_, phist::BelosMV< _ST_ > > MVT;
  typedef Belos::OperatorTraits< _ST_, phist::BelosMV< _ST_ >, st::linearOp_t > OPT;
#endif

  static void SetUpTestCase()
  {
    int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
    SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
    VTest::SetUpTestCase();
    MTest::SetUpTestCase();

    if (typeImplemented_ && !problemTooSmall_)
    {
      SUBR(linearOp_wrap_sparseMat)(&A_op,A_,&iflag_);
      ASSERT_EQ(0,iflag_);
    }
  }

  /*! Set up routine.
   */
  virtual void SetUp()
  {
    SparseMatTest::SetUp();
    VTest::SetUp();
    MTest::SetUp();

        
    if (typeImplemented_ && !problemTooSmall_)
    {
      nq_ = std::min(3*nvec_+1,(int)nglob_-4);
#ifdef HAVE_MPI
      // note: TSQR does not work if nvec>nloc (that wouldn't really be a 'tall skinny 
      // matrix' but a 'short fat and sliced matrix')
      nq_ = std::min((int)nloc_,nq_);
      int nq_local = nq_;
      iflag_ = MPI_Allreduce(&nq_local,&nq_,1,MPI_INT,MPI_MIN,mpi_comm_);
      ASSERT_EQ(0,iflag_);
#endif      
      PHISTTEST_MVEC_CREATE(&Q_,map_,nq_,&iflag_);
      ASSERT_EQ(0,iflag_);

      int v_arg[2]={_N_,_NV_};

      SUBR(mvec_put_func)(vec1_,&PHIST_TG_PREFIX(mvec123func),v_arg,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_func)(vec2_,&PHIST_TG_PREFIX(mvec321func),v_arg,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec4_,&iflag_);
      ASSERT_EQ(0,iflag_);

      // create random orthogonal Q
      SUBR(mvec_random)(Q_,&iflag_);
      ASSERT_EQ(0,iflag_);
      TYPE(sdMat_ptr) Rtmp;
      SUBR(sdMat_create)(&Rtmp,nq_,nq_,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      int rankQ = 0;
      SUBR(orthog)(NULL,Q_,NULL,Rtmp,NULL,4,&rankQ,&iflag_);
      ASSERT_GE(iflag_,0);
      SUBR(sdMat_delete)(Rtmp,&iflag_);
      ASSERT_EQ(0,iflag_);

      sigma = new _ST_[nq_];
    }

    haveMat_ = (A_ != NULL);
  }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    deleteOrthogQ();

    if( sigma != NULL)
      delete[] sigma;
    sigma = NULL;

    VTest::TearDown();
    SparseMatTest::TearDown();
    }

  void deleteOrthogQ()
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      SUBR(mvec_delete)(Q_,&iflag_); Q_=NULL;
      ASSERT_EQ(0,iflag_);
    }
  }

  static void TearDownTestCase()
  {
    MTest::TearDownTestCase();
    VTest::TearDownTestCase();
    SparseMatTest::TearDownTestCase();
  }

  // for testing the jada operator
  int nq_ = 0;
  TYPE(mvec_ptr) Q_ = NULL;
  _ST_* sigma = NULL;

  protected:

  int delete_mat(TYPE(sparseMat_ptr) A)
    {
    if (A!=NULL)
      {
      SUBR(sparseMat_delete)(A,&iflag_);
      }
    return iflag_;
    }

#ifdef PHIST_HAVE_BELOS
  int doBelosTests(TYPE(sparseMat_ptr) A)
  {
    if (typeImplemented_ && !problemTooSmall_ && haveMat_)
    {
      Teuchos::RCP<Belos::OutputManager<ST> > MyOM
        = Teuchos::rcp( new Belos::OutputManager<ST>() );
      MyOM->setVerbosity( Belos::Warnings|Belos::Debug);
      TYPE(linearOp) *op=new TYPE(linearOp);
      Teuchos::RCP<const TYPE(linearOp)> op_ptr = Teuchos::rcp(op,true);
      PHIST_ICHK_IERR(SUBR(linearOp_wrap_sparseMat)(op,A,&iflag_),iflag_);
      TYPE(linearOp) jdOp;
      // TODO setup necessary arguments for jadaOp: AX, work
      PHIST_ICHK_IERR(SUBR(jadaOp_create)(op,NULL,Q_,NULL,sigma,nq_,&jdOp,&iflag_),iflag_);
      Teuchos::RCP<const TYPE(linearOp)> jdOp_ptr=Teuchos::rcp(&jdOp,false);

      Teuchos::RCP<const phist::BelosMV< _ST_ > > V = phist::mvec_rcp< _ST_ >(vec1_,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {iflag_=-1; return iflag_;}

// note: we can't test the jadaOp in this way because it operates on a fixed number of 
// vectors and the Belos test assumes it works for any number of vectors.
//      if (Belos::TestOperatorTraits(MyOM,V,jdOp_ptr)==false) {iflag_=-2; return iflag_;}
//TODO - I'm getting an 'undefined reference' linker error here
      PHIST_ICHK_IERR(SUBR(jadaOp_delete)(&jdOp,&iflag_),iflag_);
    }
    return iflag_;
  }
#endif

  static TYPE(linearOp) A_op;

  bool haveMat_;
  };

  TYPE(linearOp) CLASSNAME::A_op;

  TEST_F(CLASSNAME, read_matrices) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_TRUE(AssertNotNull(A_));
    }
  }


#if defined(PHIST_HAVE_BELOS)

  // basic checks if the Apply function used in Belos gives the
  // same result as the standard phist sparseMat_times_mvec when
  // used with original vectors, copies or views created via the
  // Belos interface.
  TEST_F(CLASSNAME, belos_OpApply_basic)
  {
    //note: vec1==vec4 and vec2==vec3 from SetUp().
    if (typeImplemented_ && !problemTooSmall_)
    {
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1 = phist::mvec_rcp< _ST_ >(vec1_,false);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V2 = phist::mvec_rcp< _ST_ >(vec2_,false);
      
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec4_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      OPT::Apply(A_op,*V1,*V2);
      ASSERT_NEAR(1.0,MvecsEqual(vec2_,vec3_),VTest::releps());
    }
  }

  TEST_F(CLASSNAME, belos_OpApply_on_full_copy)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1 = phist::mvec_rcp< _ST_ >(vec1_,false);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V2 = phist::mvec_rcp< _ST_ >(vec2_,false);
      
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec4_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      std::vector<int> cols(nvec_);
      for (int i=0; i<nvec_; i++) cols[i]=i;
      Teuchos::RCP< phist::BelosMV< _ST_ > >  W1=MVT::CloneCopy(*V1,cols);
      Teuchos::RCP< phist::BelosMV< _ST_ > > W2=MVT::CloneCopy(*V2,cols);

      OPT::Apply(A_op,*W1,*W2);
      ASSERT_NEAR(1.0,MvecsEqual(W2->get(),vec3_),VTest::releps());
    }
  }

  TEST_F(CLASSNAME, belos_OpApply_on_colwise_copy)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1 = phist::mvec_rcp< _ST_ >(vec1_,false);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V2 = phist::mvec_rcp< _ST_ >(vec2_,false);
      
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec4_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      std::vector<int> cols(1);
      for (int i=0; i<nvec_; i++) 
      {
        cols[0]=i;
        Teuchos::RCP< phist::BelosMV< _ST_ > > W1=MVT::CloneCopy(*V1,cols);
        Teuchos::RCP< phist::BelosMV< _ST_ > > W2=MVT::CloneCopy(*V2,cols);
        OPT::Apply(A_op,*W1,*W2);
        MVT::SetBlock(*W2,cols,*V2);
      }
      ASSERT_NEAR(1.0,MvecsEqual(vec2_,vec3_),VTest::releps());
    }
  }

  TEST_F(CLASSNAME, belos_OpApply_on_full_view)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1 = phist::mvec_rcp< _ST_ >(vec1_,false);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V2 = phist::mvec_rcp< _ST_ >(vec2_,false);
      
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec4_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      std::vector<int> cols(nvec_);
      for (int i=0; i<nvec_; i++) cols[i]=i;
      Teuchos::RCP< const phist::BelosMV< _ST_ > > W1=MVT::CloneView(*V1,cols);
      Teuchos::RCP< phist::BelosMV< _ST_ > > W2=MVT::CloneViewNonConst(*V2,cols);

      OPT::Apply(A_op,*W1,*W2);
      ASSERT_NEAR(1.0,MvecsEqual(V2->get(),vec3_),VTest::releps());
      ASSERT_NEAR(1.0,MvecsEqual(W2->get(),vec3_),VTest::releps());
    }
  }

  TEST_F(CLASSNAME, belos_OpApply_on_colwise_view)
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1 = phist::mvec_rcp< _ST_ >(vec1_,false);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V2 = phist::mvec_rcp< _ST_ >(vec2_,false);
      
      SUBR(sparseMat_times_mvec)(st::one(),A_,vec4_,st::zero(),vec3_,&iflag_);
      ASSERT_EQ(0,iflag_);

      std::vector<int> cols(1);
      for (int i=0; i<nvec_; i++) 
      {
        cols[0]=i;
        Teuchos::RCP< const phist::BelosMV< _ST_ > > W1=MVT::CloneView(*V1,cols);
        Teuchos::RCP< phist::BelosMV< _ST_ > > W2=MVT::CloneViewNonConst(*V2,cols);
        OPT::Apply(A_op,*W1,*W2);
      }
      ASSERT_NEAR(1.0,MvecsEqual(vec2_,vec3_),VTest::releps());
    }
  }
  
  // e.g. in the BelosBlcoGmres solver X and Y are defined
  // to be columns/blocks j-1 and j and Y=OP*X is computed.
  // This seems to fail with GHOST on more than one MPI process,
  // so I added this extra test.
  TEST_F(CLASSNAME, belos_OpApply_on_views_of_same_vec) 
  {
    if (typeImplemented_ && !problemTooSmall_ && _NV_>1)
    {
      Teuchos::RCP<phist::BelosMV< _ST_ > > V = phist::mvec_rcp< _ST_ >(vec1_,false);
      std::vector<int> i0(1), i1(1);
      i0[0]=0; i1[0]=1;
      Teuchos::RCP<phist::BelosMV< _ST_ > > V0_copied = MVT::CloneCopy(*V,i0);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1_copied = MVT::CloneCopy(*V,i1);

      Teuchos::RCP<phist::BelosMV< _ST_ > > V0 = MVT::CloneViewNonConst(*V,i0);
      Teuchos::RCP<phist::BelosMV< _ST_ > > V1 = MVT::CloneViewNonConst(*V,i1);
      OPT::Apply(A_op,*V0_copied,*V1_copied);
      OPT::Apply(A_op,*V0,*V1);
      
      _ST_* troet;
      phist_lidx lda;

      ASSERT_REAL_EQ(1.0,MvecsEqual(V0->get(),V0_copied->get()));
      ASSERT_NEAR(1.0,MvecsEqual(V1->get(),V1_copied->get()),VTest::releps());
    }
  }

  TEST_F(CLASSNAME, belos_opTests) 
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_EQ(0,doBelosTests(A_));
    }
  }
#endif

  TEST_F(CLASSNAME,linearOp_wrap_sparseMat_apply)
  {
    if (!typeImplemented_ || problemTooSmall_)
      return;

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
  
    SUBR(sparseMat_times_mvec)(alpha,A_,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    A_op.apply(alpha,A_op.A,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec2_,vec3_));
  
  }

  TEST_F(CLASSNAME,linearOp_wrap_sparseMat_apply_fused_mvTmv)
  {
    if (!typeImplemented_ || problemTooSmall_)
      return;

    _ST_ alpha = st::prand();
    _ST_ beta = st::prand();
  
    SUBR(fused_spmv_mvTmv)(alpha,A_,vec1_,beta,vec2_,mat1_,mat2_,&iflag_);
    ASSERT_EQ(0,iflag_);

    A_op.fused_apply_mvTmv(alpha,A_op.A,vec1_,beta,vec3_,mat3_,mat4_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    ASSERT_REAL_EQ(mt::one(),MvecsEqual(vec2_,vec3_));
    ASSERT_NEAR(mt::one(),SdMatsEqual(mat1_,mat3_),std::sqrt(st::eps()));
    ASSERT_NEAR(mt::one(),SdMatsEqual(mat2_,mat4_),std::sqrt(st::eps()));
  }

  TEST_F(CLASSNAME,linearOp_wrap_sparseMat_pair_apply)
  {
    if (!typeImplemented_ || problemTooSmall_)
      return;

    TYPE(linearOp) AA_op;
    SUBR(linearOp_wrap_sparseMat_pair)(&AA_op,A_,A_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // we have v1, v2 random and v3=v2 from SetUp()

    // check that AA_op->apply is the same as A_op->apply

    _ST_ alpha=st::prand();
    _ST_ beta=st::prand();
    AA_op.apply(alpha,AA_op.A,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    A_op.apply(alpha,A_op.A,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_),VTest::releps());
        
    // clean up the operator
    AA_op.destroy(&AA_op,&iflag_);
    ASSERT_EQ(0,iflag_);
  }

  TEST_F(CLASSNAME,linearOp_wrap_sparseMat_pair_apply_shifted)
  {
    if (!typeImplemented_ || problemTooSmall_)
      return;

    TYPE(linearOp) AA_op;
    SUBR(linearOp_wrap_sparseMat_pair)(&AA_op,A_,A_,&iflag_);
    ASSERT_EQ(0,iflag_);

    // we have v1, v2 random and v3=v2 from SetUp()

    // check that  (A -1*A)x=0

    _ST_ alpha=st::prand();
    _ST_ sigma[nvec_];
    for (int i=0; i<nvec_; i++) sigma[i]=-st::one();
    AA_op.apply_shifted(alpha,AA_op.A,sigma,vec1_,st::zero(),vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),MvecEqual(vec2_,st::zero()),VTest::releps());
        
    // clean up the operator
    AA_op.destroy(&AA_op,&iflag_);
    ASSERT_EQ(0,iflag_);
  }

  TEST_F(CLASSNAME,linearOp_identity)
  {
    if (!typeImplemented_ || problemTooSmall_) return;

    TYPE(linearOp) I_op;
    SUBR(linearOp_identity)(&I_op,map_,map_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    ASSERT_TRUE(map_==I_op.range_map);
    ASSERT_TRUE(map_==I_op.domain_map);

    ST v1=st::prand();
    ST v2=st::prand();
    SUBR(mvec_put_value)(vec1_,v1,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_put_value)(vec2_,v2,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    _ST_ alpha=st::prand(), beta=st::prand();

    I_op.apply(alpha,I_op.A,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),MvecEqual(vec2_,alpha*v1+beta*v2),VTest::releps());
    
    _ST_ sigma[nvec_],asig_plusI[nvec_];
    for (int i=0; i<nvec_; i++) 
    {
      sigma[i]=-st::prand();
      asig_plusI[i]=alpha*(st::one()+sigma[i]);
    }
    
    SUBR(mvec_random)(vec1_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_random)(vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
    
    // vec3 contains alpha*sigma[j]*vec1+beta*vec2
    SUBR(mvec_vadd_mvec)(asig_plusI,vec1_,beta,vec3_,&iflag_);
    
    I_op.apply_shifted(alpha,I_op.A,sigma,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_),VTest::releps());
        
    // clean up the operator
    I_op.destroy(&I_op,&iflag_);
    ASSERT_EQ(0,iflag_);
  }

  // test wrapping AA_op=A_op*A_op (TODO: a test with two operators A!=B)
  TEST_F(CLASSNAME,linearOp_wrap_linearOp_product_apply)
  {
    if (!typeImplemented_ || problemTooSmall_)
      return;

    TYPE(linearOp) AA_op;
    SUBR(linearOp_wrap_linearOp_product)(&AA_op,&A_op,&A_op,&iflag_);
    ASSERT_EQ(0,iflag_);

    // we have v1, v2 random and v3=v2 from SetUp()

    _ST_ alpha=st::prand();
    _ST_ beta=st::prand();
    AA_op.apply(alpha,AA_op.A,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    // step-by-step to create a reference solution
    A_op.apply(alpha,A_op.A,vec1_,st::zero(),vec4_,&iflag_);
    ASSERT_EQ(0,iflag_);
    A_op.apply(st::one(),A_op.A,vec4_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);

    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_),50*VTest::releps());
        
    // clean up the operator
    SUBR(linearOp_destroy)(&AA_op,&iflag_);
    ASSERT_EQ(0,iflag_);
  }
  
    // test wrapping AAA_op=A_op*A_op*A_op and test with one Operator beeing identity 
	//(TODO: a test with two operators A!=B)
  TEST_F(CLASSNAME,linearOp_wrap_linearOp_product_triple_apply)
  {
    if (!typeImplemented_ || problemTooSmall_)
      return;

    TYPE(linearOp) AAA_op;
    SUBR(linearOp_wrap_linearOp_product_triple)(&AAA_op,&A_op,&A_op,&A_op,&iflag_);
    ASSERT_EQ(0,iflag_);

    // we have v1, v2 random and v3=v2 from SetUp()

    _ST_ alpha=st::prand();
    _ST_ beta=st::prand();
    AAA_op.apply(alpha,AAA_op.A,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
    // step-by-step to create a reference solution
    A_op.apply(alpha,A_op.A,vec1_,st::zero(),vec4_,&iflag_);
    ASSERT_EQ(0,iflag_);
	A_op.apply(st::one(),A_op.A,vec4_,st::zero(),vec5_,&iflag_);
    ASSERT_EQ(0,iflag_);
    A_op.apply(st::one(),A_op.A,vec5_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);

    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_),100*VTest::releps());
        
    // clean up the operator
    SUBR(linearOp_destroy)(&AAA_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	//test with changed identity position and linearOp_product
	TYPE(linearOp) I_op;
    SUBR(linearOp_identity)(&I_op,map_,map_,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	TYPE(linearOp) IAA_op;
    SUBR(linearOp_wrap_linearOp_product_triple)(&IAA_op,&I_op,&A_op,&A_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	TYPE(linearOp) AIA_op;
    SUBR(linearOp_wrap_linearOp_product_triple)(&AIA_op,&A_op,&I_op,&A_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	TYPE(linearOp) AAI_op;
    SUBR(linearOp_wrap_linearOp_product_triple)(&AAI_op,&A_op,&A_op,&I_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	TYPE(linearOp) AA_op;
    SUBR(linearOp_wrap_linearOp_product)(&AA_op,&A_op,&A_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
	SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec4_,&iflag_);
    ASSERT_EQ(0,iflag_);		
	SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec5_,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	IAA_op.apply(alpha,IAA_op.A,vec1_,beta,vec2_,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	AIA_op.apply(alpha,AIA_op.A,vec1_,beta,vec3_,&iflag_);
    ASSERT_EQ(0,iflag_);
	
	AAI_op.apply(alpha,AAI_op.A,vec1_,beta,vec4_,&iflag_);
    ASSERT_EQ(0,iflag_);
	
    AA_op.apply(alpha,AA_op.A,vec1_,beta,vec5_,&iflag_);
    ASSERT_EQ(0,iflag_);	
	
    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec3_),50*VTest::releps());
    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec4_),50*VTest::releps());
    ASSERT_NEAR(mt::one(),MvecsEqual(vec2_,vec5_),50*VTest::releps());
    ASSERT_NEAR(mt::one(),MvecsEqual(vec3_,vec4_),50*VTest::releps());
    ASSERT_NEAR(mt::one(),MvecsEqual(vec3_,vec5_),50*VTest::releps());
    ASSERT_NEAR(mt::one(),MvecsEqual(vec4_,vec5_),50*VTest::releps());	
	
	// clean up the operator
    SUBR(linearOp_destroy)(&IAA_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
    SUBR(linearOp_destroy)(&AIA_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
    SUBR(linearOp_destroy)(&AAI_op,&iflag_);
    ASSERT_EQ(0,iflag_);
	
    I_op.destroy(&I_op,&iflag_);
    ASSERT_EQ(0,iflag_);

    SUBR(linearOp_destroy)(&AA_op,&iflag_);
    ASSERT_EQ(0,iflag_);	
  }
  
  
