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
      read_mat("sprandn",&A1_);
      read_mat("sprandn_nodiag",&A2_);
      
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
      ASSERT_EQ(0,delete_mat(A1_));
      ASSERT_EQ(0,delete_mat(A2_));
      }
    }

TYPE(crsMat_ptr) A1_; 
TYPE(crsMat_ptr) A2_; 

protected:

int read_mat(const char* filebase,TYPE(crsMat_ptr) *ptr)
  {
  *ptr = NULL;
  char tpc = tolower(st::type_char());
  char mmfile[256],hbfile[256],binfile[256];
  sprintf(mmfile,"%c_%s%d.mm",tpc,filebase,nglob_);
#ifdef _IS_COMPLEX_
  sprintf(hbfile,"%c_%s%d.cua",tpc,filebase,nglob_);
#else  
  sprintf(hbfile,"%c_%s%d.rua",tpc,filebase,nglob_);
#endif
  sprintf(binfile,"%c_%s%d.bin",tpc,filebase,nglob_);
  
//  std::cout << "Looking for matrix \'"<<filebase<<"\'...\n";
//  std::cout << "... try \'"<<mmfile<<"\'\n";
  SUBR(crsMat_read_mm)(ptr,mmfile,&ierr_);
  if (ierr_!=_PHIST_SUCCESS_) // kernel lib can't read MatrixMarket format or file not found
    {
//    std::cout << "... try \'"<<hbfile<<"\'\n";
    SUBR(crsMat_read_hb)(ptr,hbfile,&ierr_);
    if (ierr_!=_PHIST_SUCCESS_) // kernel lib can't read Harwell-Boeing or file not found
      {
//      std::cout << "... try \'"<<binfile<<"\'\n";
      SUBR(crsMat_read_bin)(ptr,binfile,&ierr_);
      }
    }
  return ierr_;
  }

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
#ifdef PHIST_KERNEL_LIB_GHOST
      ghost_vec_t* v = (ghost_vec_t*)vec1_;
      Teuchos::RCP<const phist::GhostMV> V = phist::rcp(v,false);
#elif defined(PHIST_KERNEL_LIB_EPETRA)
      Epetra_MultiVector* v = (Epetra_MultiVector*)vec1_;
      Teuchos::RCP<const Epetra_MultiVector> V = phist::rcp(v,false);
#elif defined(PHIST_KERNEL_LIB_TPETRA)
      phist::tpetra::Traits<ST>::mvec_t* v = (phist::tpetra::Traits<ST>::mvec_t*)vec1_;
      Teuchos::RCP<const phist::tpetra::Traits<ST>::mvec_t> V = Teuchos::rcp(v,false);
#else
#endif    
      if (Belos::TestOperatorTraits(MyOM,V,op_ptr)==false) ierr_=-1;
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
