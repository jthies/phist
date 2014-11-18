// Anasazi: block krylov methods from Trilinos
void SUBR(anasazi)(TYPE(const_op_ptr) Op, 
        TYPE(mvec_ptr) vX,
        TYPE(const_mvec_ptr) vB, 
        _MT_ tol,int *num_iters, int max_blocks,
        int variant, int* nConv,
        int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#ifndef PHIST_HAVE_ANASAZI
  *ierr = -99;
#else
#include "phist_std_typedefs.hpp"  
#ifdef PHIST_KERNEL_LIB_GHOST
  typedef ghost_densemat_t MV;
  typedef phist::GhostMV AnasaziMV;
#elif defined(PHIST_KERNEL_LIB_TPETRA)
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> MV;
  typedef MV AnasaziMV;
#elif defined(PHIST_KERNEL_LIB_EPETRA)
  typedef Epetra_MultiVector MV; 
  typedef MV AnasaziMV;
#endif
  typedef st::op_t OP; // gives Sop_t, Dop_t etc.

  bool status=true;
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  int numRhs=1;// get from input vectors
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vX,&numRhs,ierr),*ierr);
  
  Teuchos::RCP<AnasaziMV> X = phist::rcp((MV*)vX, false);
  Teuchos::RCP<const AnasaziMV> B = phist::rcp((const MV*)vB, false);
  // note: our operator nas no destructor, so we should
  // actually wrap it like the ghost_densemat_t (cf. phist_rcp_helpers and phist_GhostMV),
  // but since Anasazi does nothing but apply the operator it is not necessary here.
  Teuchos::RCP<const OP> A = Teuchos::rcp((const OP*)Op, false);

////////////////////////////////////////////////////////////////////
// setup the environment, I/O, read params etc. (Teuchos package) //
////////////////////////////////////////////////////////////////////

  // create a nice wrapper for std::cout. Note that we do not allow
  // the stream to 'delete' std::cout by passing 'false' to the rcp() function
  Teuchos::RCP<Teuchos::FancyOStream> out = 
        Teuchos::rcp(new Teuchos::FancyOStream(Teuchos::rcp(&std::cout,false)));

  // customize output behavior a little:
  out->setOutputToRootOnly(0);
  out->setShowProcRank(false);

///////////////////////////////////////////////////////////////////////
// create a block GMRES "solver manager" to solve AX=B.              //
///////////////////////////////////////////////////////////////////////

  // get linear solver specific parameters from the list
  // read parameters from our XML file
  Teuchos::RCP<Teuchos::ParameterList> anasaziList 
        = Teuchos::rcp(new Teuchos::ParameterList("phist/anasazi"));
  if (variant==0 || variant==1)
  {
    //GMRES restarting
    anasaziList->set("Num Blocks",max_blocks);
  }
  if (variant==1||variant==3)
  {
    //pseudo-block, set stopping criterion
    int dq;
    if (nConv==NULL) 
    {
      dq=numRhs;
    }
    else
    {
      dq=*nConv;
    }
    anasaziList->set("Deflation Quorum",dq);
  }
  anasaziList->set("Maximum Iterations",*num_iters);
  anasaziList->set("Block Size",numRhs);
  anasaziList->set("Orthogonalization","DGKS");
  anasaziList->set("Convergence Tolerance",tol);
  anasaziList->set("Output Frequency",1);
  anasaziList->set("Show Maximum Residual Norm Only",true);
  anasaziList->set("Output Style",::Anasazi::Brief);
  int verb=::Anasazi::Errors+::Anasazi::Warnings;
#if PHIST_OUTLEV>=PHIST_INFO
  verb+=   ::Anasazi::IterationDetails
         + ::Anasazi::StatusTestDetails;
#endif
#if PHIST_OUTLEV>=PHIST_VERBOSE
  verb+=   ::Anasazi::FinalSummary
         + ::Anasazi::TimingDetails;
#endif
  anasaziList->set("Verbosity",verb);

  anasaziList->set("Output Stream",out->getOStream());

// create Anasazi problem interface
Teuchos::RCP<Anasazi::LinearProblem<ST,AnasaziMV,OP> > linearSystem
        = Teuchos::rcp(new Anasazi::LinearProblem<ST,AnasaziMV,OP>(A,X,B));

Teuchos::RCP<Anasazi::SolverManager<ST,AnasaziMV, OP> > anasazi;
if (variant==0)
  {
  anasazi = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<ST,AnasaziMV, OP>
        (linearSystem, anasaziList));
  }
else
  {
  PHIST_OUT(PHIST_ERROR,"anasazi variant %d is not supported",variant);
  *ierr=-99;
  }

#if PHIST_OUTLEV>PHIST_DEBUG
  PHIST_DEB("valid Anasazi parameters:\n");
  std::cerr << *anasazi->getValidParameters();
  PHIST_DEB("current GMRES parameters:\n");
  std::cerr << *anasazi->getCurrentParameters();
#endif

///////////////////////////////////////////////////////////////////////
// solve the system!                                                 //
///////////////////////////////////////////////////////////////////////
try {
    linearSystem->setProblem();
    Anasazi::ReturnType result=anasazi->solve();
    *num_iters = anasazi->getNumIters();
    *out << "Anasazi returned '"<<Anasazi::convertReturnTypeToString(result)<<"'"<<std::endl;
    if (result!=Anasazi::Converged) *ierr=1;
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,status);

  if (!status) *ierr=PHIST_CAUGHT_EXCEPTION; 
  return;
#endif /* PHIST_HAVE_BELOS */
  }// end of anasazi
