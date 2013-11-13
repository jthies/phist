// block GMRES - solve a linear system with multiple RHS
// using the BlockGMRES method implemented in Belos.
// AX=B is solved for X. At most max_blocks blocks of vectors are
// generated. A block consists of as many vectors as there
// are columns in X and B. If max_blocks is reached, the  
// method is restarted, until num_iters is reached or     
// ||r||/||r0||<tol is achieved. *num_iters is overwritten
// by the actual number of iterations performed.
void SUBR(bgmres)(TYPE(const_op_ptr) Op, 
        TYPE(mvec_ptr) vX,
        TYPE(const_mvec_ptr) vB, 
        _MT_ tol,int *num_iters, int max_blocks,
        int variant, int* nConv,
        int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
#ifdef PHIST_KERNEL_LIB_GHOST
  typedef ghost_vec_t MV;
  typedef phist::GhostMV BelosMV;
#elif defined(PHIST_KERNEL_LIB_TPETRA)
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> MV;
  typedef MV BelosMV;
#elif defined(PHIST_KERNEL_LIB_EPETRA)
  typedef Epetra_MultiVector MV; 
  typedef MV BelosMV;
#endif
  typedef st::op_t OP; // gives Sop_t, Dop_t etc.

  bool status=true;
  *ierr=0;
  
  int numRhs=1;// get from input vectors
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vX,&numRhs,ierr),*ierr);
  
  Teuchos::RCP<BelosMV> X = phist::rcp((MV*)vX, false);
  Teuchos::RCP<const BelosMV> B = phist::rcp((const MV*)vB, false);
  // note: our operator nas no destructor, so we should
  // actually wrap it like the ghost_vec_t (cf. phist_rcp_helpers and phist_GhostMV),
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

  belosList->set("Maximum Iterations",*num_iters);
  belosList->set("Block Size",numRhs);
  belosList->set("Num Blocks",max_blocks);
  belosList->set("Orthogonalization","DGKS");
  belosList->set("Convergence Tolerance",tol);
  belosList->set("Output Frequency",1);
  belosList->set("Show Maximum Residual Norm Only",true);
  belosList->set("Output Style",::Belos::Brief);
  belosList->set("Verbosity",::Belos::Errors+::Belos::Warnings
                         +::Belos::IterationDetails
                         +::Belos::StatusTestDetails
                         +::Belos::FinalSummary
                         +::Belos::TimingDetails);

  belosList->set("Output Stream",out->getOStream());

// create Belos problem interface
Teuchos::RCP<Belos::LinearProblem<ST,BelosMV,OP> > linearSystem
        = Teuchos::rcp(new Belos::LinearProblem<ST,BelosMV,OP>(A,X,B));

Teuchos::RCP<Belos::SolverManager<ST,BelosMV, OP> > gmres;
if (variant==0)
  {
  gmres = Teuchos::rcp(new Belos::BlockGmresSolMgr<ST,BelosMV, OP>
        (linearSystem, belosList));
  }
else if (variant==1)
  {
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
  gmres = Teuchos::rcp(new Belos::PseudoBlockGmresSolMgr<ST,BelosMV, OP>
        (linearSystem, belosList));
  }
else
  {
  PHIST_OUT(PHIST_ERROR,"gmres variant %d is not supported",variant);
  *ierr=-99;
  }

#if PHIST_OUTLEV>PHIST_DEBUG
  PHIST_DEB("valid GMRES parameters:");
  std::cerr << *gmres->getValidParameters();
  PHIST_DEB("current GMRES parameters:");
  std::cerr << *gmres->getCurrentParameters();
#endif

///////////////////////////////////////////////////////////////////////
// solve the system!                                                 //
///////////////////////////////////////////////////////////////////////
try {
    linearSystem->setProblem();
    Belos::ReturnType result=gmres->solve();
    *num_iters = gmres->getNumIters();
    *out << "Belos returned '"<<Belos::convertReturnTypeToString(result)<<"'"<<std::endl;
    if (result!=Belos::Converged) *ierr=1;
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,status);

  if (!status) *ierr=PHIST_CAUGHT_EXCEPTION; 
  return;
  }// end of bgmres


  
