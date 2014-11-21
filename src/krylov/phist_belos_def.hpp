#if defined PHIST_KERNEL_LIB_FORTRAN
// use the correct rcp function, this is a bit adhoc
// and should be implemented nicer
#define PHIST_RCP phist::PREFIX(rcp)
#else
#define PHIST_RCP phist::rcp
#endif

// Belos: block krylov methods from Trilinos
extern "C" void SUBR(belos)(TYPE(const_op_ptr) Op, 
        TYPE(mvec_ptr) vX,
        TYPE(const_mvec_ptr) vB, 
        _MT_ tol,int *num_iters, int max_blocks,
        int variant, int* nConv,
        int* ierr)
  {
#ifndef PHIST_HAVE_BELOS
  ENTER_FCN(__FUNCTION__);
  *ierr = -99;
#else
#include "phist_std_typedefs.hpp"  
#ifdef PHIST_KERNEL_LIB_GHOST
  typedef ghost_densemat_t MV;
  typedef phist::GhostMV BelosMV;
#elif defined(PHIST_KERNEL_LIB_TPETRA)
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> MV;
  typedef MV BelosMV;
#elif defined(PHIST_KERNEL_LIB_EPETRA)
  typedef Epetra_MultiVector MV; 
  typedef MV BelosMV;
#else
typedef st::mvec_t MV;
typedef phist::MultiVector< _ST_ > BelosMV;
#endif
  typedef st::op_t OP; // gives Sop_t, Dop_t etc.

  bool status=true;
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  int numRhs=1;// get from input vectors
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vX,&numRhs,ierr),*ierr);
  
  Teuchos::RCP<BelosMV> X = PHIST_RCP((MV*)vX, false);
  Teuchos::RCP<const BelosMV> B = PHIST_RCP((const MV*)vB, false);
  // note: our operator nas no destructor, so we should
  // actually wrap it like the ghost_densemat_t (cf. phist_rcp_helpers and phist_GhostMV),
  // but since Belos does nothing but apply the operator it is not necessary here.
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
  Teuchos::RCP<Teuchos::ParameterList> belosList 
        = Teuchos::rcp(new Teuchos::ParameterList("phist/belos"));
  if (variant==0 || variant==1)
  {
    //GMRES restarting
    belosList->set("Num Blocks",max_blocks);
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
    belosList->set("Deflation Quorum",dq);
  }
  belosList->set("Maximum Iterations",*num_iters);
  belosList->set("Block Size",numRhs);
  belosList->set("Orthogonalization","DGKS");
  belosList->set("Convergence Tolerance",tol);
  belosList->set("Output Frequency",1);
  belosList->set("Show Maximum Residual Norm Only",true);
  belosList->set("Output Style",::Belos::Brief);
  int verb=::Belos::Errors+::Belos::Warnings;
#if PHIST_OUTLEV>=PHIST_INFO
  verb+=   ::Belos::IterationDetails
         + ::Belos::StatusTestDetails;
#endif
#if PHIST_OUTLEV>=PHIST_VERBOSE
  verb+=   ::Belos::FinalSummary
         + ::Belos::TimingDetails;
#endif
  belosList->set("Verbosity",verb);

  belosList->set("Output Stream",out->getOStream());

// create Belos problem interface
Teuchos::RCP<Belos::LinearProblem<ST,BelosMV,OP> > linearSystem
        = Teuchos::rcp(new Belos::LinearProblem<ST,BelosMV,OP>(A,X,B));

Teuchos::RCP<Belos::SolverManager<ST,BelosMV, OP> > belos;
if (variant==0)
  {
  belos = Teuchos::rcp(new Belos::BlockGmresSolMgr<ST,BelosMV, OP>
        (linearSystem, belosList));
  }
else if (variant==1)
  {
  belos = Teuchos::rcp(new Belos::PseudoBlockGmresSolMgr<ST,BelosMV, OP>
        (linearSystem, belosList));
  }
else if (variant==2)
  {
  belos = Teuchos::rcp(new Belos::BlockCGSolMgr<ST,BelosMV, OP>
        (linearSystem, belosList));
  }
else if (variant==3)
  {
  belos = Teuchos::rcp(new Belos::PseudoBlockCGSolMgr<ST,BelosMV, OP>
        (linearSystem, belosList));
  }
else
  {
  PHIST_OUT(PHIST_ERROR,"belos variant %d is not supported",variant);
  *ierr=-99;
  }

#if PHIST_OUTLEV>PHIST_DEBUG
  PHIST_DEB("valid Belos parameters:\n");
  std::cerr << *belos->getValidParameters();
  PHIST_DEB("current GMRES parameters:\n");
  std::cerr << *belos->getCurrentParameters();
#endif

///////////////////////////////////////////////////////////////////////
// solve the system!                                                 //
///////////////////////////////////////////////////////////////////////
try {
    linearSystem->setProblem();
    Belos::ReturnType result=belos->solve();
    *num_iters = belos->getNumIters();
    *out << "Belos returned '"<<Belos::convertReturnTypeToString(result)<<"'"<<std::endl;
    if (result!=Belos::Converged) *ierr=1;
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,status);

  if (!status) *ierr=PHIST_CAUGHT_EXCEPTION; 
  return;
#endif /* PHIST_HAVE_BELOS */
  }// end of belos

#undef PHIST_RCP
