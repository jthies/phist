#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,1> 
  {

public:

  typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
  typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,1> VTest;

  static void SetUpTestCase()
  {
    SparseMatTest::SetUpTestCase();
    VTest::SetUpTestCase();
  }

  /*! Set up routine.
   */
  virtual void SetUp()
    {
    SparseMatTest::SetUp();
    VTest::SetUp();

    createOrthogQ();
    
    if (typeImplemented_ && !problemTooSmall_)
      {
      sigma = new _ST_[nq_];
      }

    haveMat_ = (A_ != NULL);
    }

  void createOrthogQ()
  {
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
    }
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
      TYPE(op) *op=new TYPE(op);
      Teuchos::RCP<const TYPE(op)> op_ptr = Teuchos::rcp(op,true);
      PHIST_ICHK_IERR(SUBR(op_wrap_sparseMat)(op,A,&iflag_),iflag_);
      TYPE(op) jdOp;
      // TODO setup necessary arguments for jadaOp: AX, work
      PHIST_ICHK_IERR(SUBR(jadaOp_create)(op,NULL,Q_,NULL,sigma,nq_,&jdOp,&iflag_),iflag_);
      Teuchos::RCP<const TYPE(op)> jdOp_ptr=Teuchos::rcp(&jdOp,false);
#ifdef PHIST_KERNEL_LIB_GHOST
      ghost_densemat_t* v = (ghost_densemat_t*)vec1_;
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

  bool haveMat_;
  };

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
//TODO test jada operator, test operator and compare with raw CrsMat
