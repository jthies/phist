#ifndef CLASSNAME
#error "file not included correctly."
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
    
    if (typeImplemented_)
      {
      nq_ = std::min(3*nvec_+1,nglob_-4);
#ifdef HAVE_MPI
      // note: TSQR does not work if nvec>nloc (that wouldn't really be a 'tall skinny 
      // matrix' but a 'short fat and sliced matrix')
      nq_ = std::min(nloc_,nq_);
      int nq_local = nq_;
      ierr_ = MPI_Allreduce(&nq_local,&nq_,1,MPI_INT,MPI_MIN,mpi_comm_);
      ASSERT_EQ(0,ierr_);
#endif      
      SUBR(mvec_create)(&Q_,map_,nq_,&ierr_);
      ASSERT_EQ(0,ierr_);
      sigma = new _ST_[nq_];

      // create random orthogonal Q
      SUBR(mvec_random)(Q_,&ierr_);
      ASSERT_EQ(0,ierr_);
      TYPE(sdMat_ptr) Rtmp;
      SUBR(sdMat_create)(&Rtmp,nq_,nq_,comm_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_QR)(Q_,Rtmp,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(sdMat_delete)(Rtmp,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(read_mat)("sprandn",nglob_,&A1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(read_mat)("sprandn_nodiag",nglob_,&A2_,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      if (A1_==NULL || A2_==NULL)
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
    SUBR(mvec_delete)(Q_,&ierr_);
    ASSERT_EQ(0,ierr_);
    delete[] sigma;
      ASSERT_EQ(0,delete_mat(A1_));
      ASSERT_EQ(0,delete_mat(A2_));
      }
    }

TYPE(crsMat_ptr) A1_; 
TYPE(crsMat_ptr) A2_; 

// for testing the jada operator
int nq_;
TYPE(mvec_ptr) Q_;
_ST_* sigma;

protected:

int delete_mat(TYPE(crsMat_ptr) A)
  {
  if (A!=NULL)
    {
    SUBR(crsMat_delete)(A,&ierr_);
    }
  return ierr_;
  }

int doBelosTests(TYPE(crsMat_ptr) A)
  {
    if (typeImplemented_ && haveMats_)
      {
      Teuchos::RCP<Belos::OutputManager<ST> > MyOM
        = Teuchos::rcp( new Belos::OutputManager<ST>() );
      MyOM->setVerbosity( Belos::Warnings|Belos::Debug);
      TYPE(op) *op=new TYPE(op);
      Teuchos::RCP<const TYPE(op)> op_ptr = Teuchos::rcp(op,true);
      PHIST_ICHK_IERR(SUBR(op_wrap_crsMat)(op,A,&ierr_),ierr_);
      TYPE(op) jdOp;
      // TODO setup necessary arguments for jadaOp: AX, work
      PHIST_ICHK_IERR(SUBR(jadaOp_create)(op,NULL,Q_,NULL,sigma,_NV_,&jdOp,&ierr_),ierr_);
      Teuchos::RCP<const TYPE(op)> jdOp_ptr=Teuchos::rcp(&jdOp,false);
#ifdef PHIST_KERNEL_LIB_GHOST
      ghost_densemat_t* v = (ghost_densemat_t*)vec1_;
      Teuchos::RCP<const phist::GhostMV> V = phist::rcp(v,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {ierr_=-1; return ierr_;}
#elif defined(PHIST_KERNEL_LIB_EPETRA)
      Epetra_MultiVector* v = (Epetra_MultiVector*)vec1_;
      Teuchos::RCP<const Epetra_MultiVector> V = phist::rcp(v,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {ierr_=-1; return ierr_;}
#elif defined(PHIST_KERNEL_LIB_TPETRA)
      phist::tpetra::Traits<ST>::mvec_t* v = (phist::tpetra::Traits<ST>::mvec_t*)vec1_;
      Teuchos::RCP<const phist::tpetra::Traits<ST>::mvec_t> V = Teuchos::rcp(v,false);
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) {ierr_=-1; return ierr_;}
#else
#warning belos test case not implemented for this kernel lib (OpTest_def_hpp)
#endif
// note: we can't test the jadaOp in this way because it operates on a fixed number of 
// vectors and the Belos test assumes it works for any number of vectors.
//      if (Belos::TestOperatorTraits(MyOM,V,jdOp_ptr)==false) {ierr_=-2; return ierr_;}
//TODO - I'm getting an 'undefined reference' linker error here
      PHIST_ICHK_IERR(SUBR(jadaOp_delete)(&jdOp,&ierr_),ierr_);
      }
  return ierr_;
  }

bool haveMats_;
};

  TEST_F(CLASSNAME, read_matrices) 
    {
    if (typeImplemented_)
      {
      ASSERT_TRUE(AssertNotNull(A1_));
      ASSERT_TRUE(AssertNotNull(A2_));
      }
    }


  TEST_F(CLASSNAME, belos_opTests) 
    {
    if (typeImplemented_)
      {
      ASSERT_EQ(0,doBelosTests(A1_));
      ASSERT_EQ(0,doBelosTests(A2_));
      }
    }

//TODO test jada operator, test operator and compare with raw CrsMat
