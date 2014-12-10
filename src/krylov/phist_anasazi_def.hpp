#ifdef PHIST_KERNEL_LIB_BUILTIN
// this should be done for wrapping the actual phist void-pointers
// as typed multivectors, this macro will expand to S/D/C/Zrcp because
// overloading is not possible on return type only.
#define PHIST_rcp phist::PREFIX(rcp)
#else
// For C++ objects (E/Tpetra libs) this results in Teuchos::rcp,
// for GHOST we have our own implementation of the function
#define PHIST_rcp phist::rcp
#endif
// Anasazi: block krylov methods from Trilinos
void SUBR(anasazi)(      TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) Ainv_op, 
                         TYPE(const_op_ptr) B_op,
                         TYPE(const_mvec_ptr) v0,  eigSort_t which,
                         _MT_ tol,                 int *nEig,
                         int* nIter,               int blockDim,
                         int numBlocks,
                         bool symmetric,
                         TYPE(mvec_ptr) vX,         _ST_* eigs,
                         int* ierr)
  {
  PHIST_ENTER_FCN(__FUNCTION__);
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
#else
  typedef st::mvec_t MV;
  typedef phist::MultiVector< _ST_ > AnasaziMV;
#endif
  typedef st::op_t OP; // gives Sop_t, Dop_t etc.
  typedef Anasazi::MultiVecTraits<ST,AnasaziMV> MVT;

  bool status=true;
  *ierr=0;
  
  int variant=0;// only block Krylov-Schur implemented
  
  Teuchos::RCP<AnasaziMV> X = PHIST_rcp((MV*)vX, false);
  Teuchos::RCP<AnasaziMV> X0;
  if (v0!=NULL)
  {
    X0 = PHIST_rcp((MV*)v0, false);
  }
  else
  {
    X0=MVT::Clone(*X,blockDim);
  }

  // note: our operator nas no destructor, so we should
  // actually wrap it like the ghost_densemat_t (cf. phist_rcp_helpers and phist_GhostMV),
  // but since Anasazi does nothing but apply the operator it is not necessary here.
  Teuchos::RCP<const OP> A = Teuchos::rcp((const OP*)A_op, false);
  Teuchos::RCP<const OP> B=Teuchos::null;
  if (B_op!=NULL)
  {
     B = Teuchos::rcp((const OP*)B_op, false);
  }

  // this can be used to provide e.g. the operation (A-sigma*B)^{-1}
  Teuchos::RCP<const OP> Op=Teuchos::null;

  if (Ainv_op!=NULL)
  {
     Op = Teuchos::rcp((const OP*)Ainv_op, false);
  }

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

  anasaziList->set("Which",std::string(eigSort2str(which)));
  anasaziList->set("Block Size",blockDim);
  anasaziList->set("Num Blocks",numBlocks);
  anasaziList->set("Maximum Restarts",(int)(*nIter/numBlocks));
  anasaziList->set("Orthogonalization","SVQB");
  anasaziList->set("Convergence Tolerance",tol);
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
  
  Teuchos::RCP<Anasazi::Eigenproblem<ST,AnasaziMV,OP> > eigenProblem
        = Teuchos::rcp(new Anasazi::BasicEigenproblem<ST,AnasaziMV,OP>());

eigenProblem->setA(A);
if (B!=Teuchos::null)
{
  eigenProblem->setM(B);
}
if (Op!=Teuchos::null)
{
  eigenProblem->setOperator(Op);
}
eigenProblem->setHermitian(symmetric);        
eigenProblem->setNEV(*nEig);
eigenProblem->setInitVec(X0);

    PHIST_CHK_IERR(*ierr=eigenProblem->setProblem()?0:-1,*ierr);

Teuchos::RCP<Anasazi::SolverManager<ST,AnasaziMV, OP> > anasazi;
if (variant==0)
  {
  anasazi = Teuchos::rcp(new Anasazi::BlockKrylovSchurSolMgr<ST,AnasaziMV, OP>
        (eigenProblem, *anasaziList));
  }
else
  {
  PHIST_OUT(PHIST_ERROR,"anasazi variant %d is not supported",variant);
  *ierr=-99;
  }

///////////////////////////////////////////////////////////////////////
// solve the system!                                                 //
///////////////////////////////////////////////////////////////////////
try {
    Anasazi::ReturnType result=anasazi->solve();
    *nIter = anasazi->getNumIters();
    *out << "Anasazi returned '"
        << (result==Anasazi::Converged?  "Converged":
          (result==Anasazi::Unconverged?"Unconverged":
                               "<unknown return flag>"))
          <<"'"<<std::endl;
    if (result!=Anasazi::Converged)
    { 
      PHIST_SOUT(PHIST_WARNING,"if not all requested eigenvalues\n"
                               "converged, the returned eigenvectors\n"
                               "are incorrect, probably due to incorrect\n"
                               "sorting or a missing function call (issue #79)\n");
      *ierr=1;
    }
  } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,*out,status);

  if (!status)
  {
    *ierr=PHIST_CAUGHT_EXCEPTION;
    return;
  }
  
  const Anasazi::Eigensolution<ST,AnasaziMV>& soln=eigenProblem->getSolution();
  *nEig=soln.numVecs;
#ifndef IS_COMPLEX
  if (!symmetric)
  {
    PHIST_SOUT(PHIST_WARNING,"%s, returning only real parts of eigenvalues\n"
                             "and the eigenvectors are stored in compressed form\n"
                             "(file %s, line %d)\n",
                             __FUNCTION__,__FILE__,__LINE__);
  }
#endif
  for (int i=0;i<*nEig; i++)
  {
#ifdef IS_COMPLEX
  eigs[i]=std::complex<MT>(soln.Evals[i].realpart,soln.Evals[i].imagpart);
#else
  eigs[i]=soln.Evals[i].realpart;
#endif
  }

  if (*nEig>0)
  {
    AnasaziMV* evecs = soln.Evecs.getRawPtr();
    MV* _evecs=(MV*)evecs;
#if defined(PHIST_KERNEL_LIB_GHOST)||defined(PHIST_KERNEL_LIB_FORTRAN)
    _evecs=evecs->get();
#endif
    PHIST_CHK_IERR(SUBR(mvec_get_block)(_evecs,vX,0,*nEig-1, ierr),*ierr);
  }
  PHIST_SOUT(PHIST_VERBOSE,"returning nEig=%d, nIter=%d\n",*nEig,*nIter);
  return;
#endif /* PHIST_HAVE_ANASAZI */
  }// end of anasazi
#undef PHIST_rcp
