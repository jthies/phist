#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,3>, 
                 public virtual KernelTestWithSdMats<_ST_,_NV_,_NV_,0,4>
  {

public:

  typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,3> VTest;
  typedef KernelTestWithSdMats<_ST_,_NV_,_NV_,0,4> MTest;

  static void SetUpTestCase()
  {
    int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
    SparseMatTest::SetUpTestCase(sparseMatCreateFlag);
    VTest::SetUpTestCase();
    MTest::SetUpTestCase();

    SUBR(linearOp_wrap_sparseMat)(&A_op,A_,&iflag_);
    ASSERT_EQ(0,iflag_);
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

      SUBR(mvec_random)(vec1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_random)(vec2_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),vec2_,st::zero(),vec3_,&iflag_);
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
#ifdef PHIST_KERNEL_LIB_GHOST
      ghost_densemat* v = (ghost_densemat*)vec1_;
      Teuchos::RCP<const phist::GhostMV> V = phist::rcp(v,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {iflag_=-1; return iflag_;}
#elif defined(PHIST_KERNEL_LIB_EPETRA)
      Epetra_MultiVector* v = (Epetra_MultiVector*)vec1_;
      Teuchos::RCP<const Epetra_MultiVector> V = phist::rcp(v,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {iflag_=-1; return iflag_;}
#elif defined(PHIST_KERNEL_LIB_TPETRA)
      phist::tpetra::Traits<ST>::mvec_t* v = (phist::tpetra::Traits<ST>::mvec_t*)vec1_;
      Teuchos::RCP<const phist::tpetra::Traits<ST>::mvec_t> V = Teuchos::rcp(v,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {iflag_=-1; return iflag_;}
#else
#warning belos test case not implemented for this kernel lib (OpTest_def_hpp)
#endif
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
#ifdef PHIST_KERNEL_LIB_GHOST
  TEST_F(CLASSNAME, DISABLED_belos_opTests) 
#else
  TEST_F(CLASSNAME, belos_opTests) 
#endif
  {
    if (typeImplemented_ && !problemTooSmall_)
    {
      ASSERT_EQ(0,doBelosTests(A_));
    }
  }
#endif

  TEST_F(CLASSNAME,linearOp_wrap_sparseMat_apply)
  {
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
