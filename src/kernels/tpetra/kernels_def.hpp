#include "Teuchos_StandardCatchMacros.hpp"
#include "../cpp_macros.h"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_MatrixIO.hpp"

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
_SUBROUTINE_(read_crsMat_mm)(void** A, char* filename,int* ierr)
  {
  // TODO - doesn't seem to be available in Tpetra
  *ierr=-99;
  }

//! read a matrix from a Ghost CRS (binary) file.
_SUBROUTINE_(read_crsMat_bin)(void** A, char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
_SUBROUTINE_(read_crsMat_hb)(void** vA, char* filename,int* ierr)
  {
  *ierr=0;
  std::string fname(filename);
  Teuchos::RCP<const comm_t> comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::ParameterList nodeParams;
  
  Teuchos::RCP<Traits<_ST_>::crsMat_t> A;
  Teuchos::RCP<node_t> node = Teuchos::rcp(new node_t(nodeParams));
  Teuchos::RCP<map_t> rowMap=Teuchos::null; // assume linear distribution for now
  Teuchos::RCP<Teuchos::ParameterList> params=Teuchos::null;
  //_TRY_CATCH_(Tpetra::Utils::readHBMatrix(fname,comm,node,A,rowMap, params),*ierr);
  _TRY_CATCH_(Tpetra::Utils::readHBMatrix(fname,comm,node,A),*ierr);
  *vA = (void*)(A.get());
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
_SUBROUTINE_(get_row_map)(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (void*)(A->getRowMap().get());
  }

//! get column distribution of a matrix
_SUBROUTINE_(get_col_map)(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (void*)(A->getColMap().get());
  }

//! get the map for vectors x in y=A*x
_SUBROUTINE_(get_domain_map)(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (void*)(A->getDomainMap().get());
  }

//! get the map for vectors y in y=A*x
_SUBROUTINE_(get_range_map)(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (void*)(A->getRangeMap().get());
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
_SUBROUTINE_(create_mvec)(void* vmap, int nvec, void** vV, _ST_** val, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(map_t, map, vmap, *ierr);
  Teuchos::RCP<map_t> map_ptr = Teuchos::rcp(map,false);
  Traits<_ST_>::mvec_t* result = new Traits<_ST_>::mvec_t(map_ptr,nvec);
  vV=(void**)(&result);
  Teuchos::ArrayRCP<_ST_> val_ptr = result->get1dViewNonConst();
  *val = val_ptr.getRawPtr();
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
_SUBROUTINE_(create_sdMat)(int nrows, int ncols, void** vM, _ST_** val, int* ierr)
  {
  *ierr=0;
  // create local map
  Teuchos::RCP<Teuchos::SerialComm<int> > scomm = Teuchos::rcp
        (new Teuchos::SerialComm<int>());
  //TODO - add node arg
  Teuchos::RCP<map_t> localMap =
        Teuchos::rcp(new map_t(nrows, 0, scomm, Tpetra::LocallyReplicated));
  Traits<_ST_>::sdMat_t* result = new Traits<_ST_>::mvec_t(localMap,ncols);
  *vM=(void*)(result);
  Teuchos::ArrayRCP<_ST_> valptr = result->get1dViewNonConst();
  *val = valptr.getRawPtr();
  }

//@}

//! \name destructors

//@{

//!
_SUBROUTINE_(delete_crsMat)(void* vA, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::crsMat_t,A,vA,*ierr);
  delete [] A;
  }

//!
_SUBROUTINE_(delete_mvec)(void* vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  delete [] V;
  }

//!
_SUBROUTINE_(delete_sdMat)(void* vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,M,vM,*ierr);
  delete [] M;
  }

//@}

//! \name Numerical functions
//!@{

//! y=alpha*A*x+beta*y.
_SUBROUTINE_(crsMat_X_mvec)(_ST_ alpha, void* vA, void* vx, _ST_ beta, void* vy, int* ierr)
  {
  *ierr=0;
  
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t,A,vA,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,x,vx,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,y,vy,*ierr);

  Traits<_ST_>::crsMVM_t spMVM(Teuchos::rcp(A,false));
  _TRY_CATCH_(spMVM.apply(*x,*y,Teuchos::NO_TRANS,alpha,beta),*ierr);

  }

//! dot product of vectors v_i and w_i, i=1..numvecs
_SUBROUTINE_(mvec_dot_mvec)(void* vv, void* vw, _ST_* s, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,v,vv,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,w,vw,*ierr);
  Teuchos::ArrayView<_ST_> dots(s,v->getNumVectors());
  _TRY_CATCH_(v->dot(*w,dots),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
_SUBROUTINE_(mvecT_X_mvec)(void* vV, void* vW, void* vC, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,C,vC,*ierr);
  _ST_ alpha = Traits<_ST_>::one();
  _ST_ beta = Traits<_ST_>::zero();
  _TRY_CATCH_(C->multiply(Teuchos::TRANS, Teuchos::NO_TRANS,alpha,*V,*W,beta),*ierr);
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
_SUBROUTINE_(mvec_QR)(void* V, void* Q, void* R, int* ierr)
  {
  // TODO - check the status of TSQR in Trilinos
  *ierr=-99;
  }

//!@}
