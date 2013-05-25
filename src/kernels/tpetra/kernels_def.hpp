#include "phist_macros.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "../cpp_macros.h"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_MatrixIO.hpp"

extern "C" {

// we implement all the four types
void _SUBR_(type_avail)(int* ierr)
  {
  *ierr=0;
  }

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void _SUBR_(crsMat_read_mm)(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  // TODO - doesn't seem to be available in Tpetra
  *ierr=-99;
  }

//! read a matrix from a Ghost CRS (binary) file.
void _SUBR_(crsMat_read_bin)(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void _SUBR_(crsMat_read_hb)(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
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
  *vA = (_TYPE_(crsMat_ptr))(A.get());
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void _SUBR_(crsMat_get_row_map)(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* 
ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getRowMap().get());
  }

//! get column distribution of a matrix
void _SUBR_(crsMat_get_col_map)(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getColMap().get());
  }

//! get the map for vectors x in y=A*x
void _SUBR_(crsMat_get_domain_map)(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getDomainMap().get());
  }

//! get the map for vectors y in y=A*x
void _SUBR_(crsMat_get_range_map)(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getRangeMap().get());
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void _SUBR_(mvec_create)(_TYPE_(mvec_ptr)* vV, const_map_ptr_t vmap, int nvec, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t, map, vmap, *ierr);
  Teuchos::RCP<const map_t> map_ptr = Teuchos::rcp(map,false);
  Traits<_ST_>::mvec_t* V = new Traits<_ST_>::mvec_t(map_ptr,nvec);
  *vV=(_TYPE_(mvec_ptr))(V);
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void _SUBR_(sdMat_create)(_TYPE_(sdMat_ptr)* vM, int nrows, int ncols, int* ierr)
  {
  *ierr=0;
  // create local map
  Teuchos::RCP<Teuchos::SerialComm<int> > scomm = Teuchos::rcp
        (new Teuchos::SerialComm<int>());
  //TODO - add node arg
  Teuchos::RCP<map_t> localMap =
        Teuchos::rcp(new map_t(nrows, 0, scomm, Tpetra::LocallyReplicated));
  Traits<_ST_>::sdMat_t* M = new Traits<_ST_>::mvec_t(localMap,ncols);
  *vM=(_TYPE_(sdMat_ptr))(M);
  }

//@}

//! retrieve local length of the vectors in V
void _SUBR_(mvec_my_length)(_TYPE_(const_mvec_ptr) vV, int* len, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  *len = V->getLocalLength();
  }

//! retrieve number of vectors/columns in V
void _SUBR_(mvec_num_vectors)(_TYPE_(const_mvec_ptr) vV, int* nvec, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  *nvec = V->getNumVectors();
  }


//! extract view from multi-vector
void _SUBR_(mvec_extract_view)(_TYPE_(mvec_ptr) vV, _ST_** val, int* lda, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  Teuchos::ArrayRCP<_ST_> val_ptr = V->get1dViewNonConst();
  *val = val_ptr.getRawPtr();
  *lda = V->getLocalLength();
  }

//! extract view from serial dense matrix
void _SUBR_(sdMat_extract_view)(_TYPE_(sdMat_ptr) vM, _ST_** val, int* lda, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,M,vM,*ierr);
  Teuchos::ArrayRCP<_ST_> valptr = M->get1dViewNonConst();
  *val = valptr.getRawPtr();
  *lda=0;
  *ierr=-99;
  }


//! \name destructors

//@{

//!
void _SUBR_(crsMat_delete)(_TYPE_(crsMat_ptr) vA, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::crsMat_t,A,vA,*ierr);
  delete A;
  }

//!
void _SUBR_(mvec_delete)(_TYPE_(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  delete V;
  }

//!
void _SUBR_(sdMat_delete)(_TYPE_(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,M,vM,*ierr);
  delete M;
  }

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
void _SUBR_(mvec_put_value)(_TYPE_(mvec_ptr) vV, _ST_ value, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _TRY_CATCH_(V->putScalar(value),*ierr);
  }

//! put random numbers into all elements of a multi-vector
void _SUBR_(mvec_random)(_TYPE_(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _TRY_CATCH_(V->randomize(),*ierr);
  }

//! put random numbers into all elements of a serial dense matrix
void _SUBR_(sdMat_random)(_TYPE_(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,M,vM,*ierr);
  _TRY_CATCH_(M->randomize(),*ierr);
  }


//! y=alpha*x+beta*y
void _SUBR_(mvec_add_mvec)(_ST_ alpha, _TYPE_(const_mvec_ptr) vX,
                            _ST_ beta,  _TYPE_(mvec_ptr)       vY, 
                            int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,X,vX,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,Y,vX,*ierr);
  _TRY_CATCH_(Y->update(alpha,*X,beta),*ierr);
  }


//! y=alpha*A*x+beta*y.
void _SUBR_(crsMat_times_mvec)(_ST_ alpha, _TYPE_(const_crsMat_ptr) vA, 
                                        _TYPE_(const_mvec_ptr) vx, 
                                        _ST_ beta, _TYPE_(mvec_ptr) vy, 
                                        int* ierr)
  {
  *ierr=0;
  
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t,A,vA,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,x,vx,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,y,vy,*ierr);

  Traits<_ST_>::crsMVM_t spMVM(Teuchos::rcp(A,false));
  _TRY_CATCH_(spMVM.apply(*x,*y,Teuchos::NO_TRANS,alpha,beta),*ierr);

  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void _SUBR_(mvec_dot_mvec)(_TYPE_(const_mvec_ptr) vv, _TYPE_(const_mvec_ptr) vw, _ST_* s, 
int* 
ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,v,vv,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,w,vw,*ierr);
  Teuchos::ArrayView<_ST_> dots(s,v->getNumVectors());
  _TRY_CATCH_(v->dot(*w,dots),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void _SUBR_(mvecT_times_mvec)(_ST_ alpha, _TYPE_(const_mvec_ptr) vV, 
                           _TYPE_(const_mvec_ptr) vW, _ST_ beta, 
                           _TYPE_(sdMat_ptr) vC, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,C,vC,*ierr);
  _TRY_CATCH_(C->multiply(Teuchos::TRANS, Teuchos::NO_TRANS,alpha,*V,*W,beta),*ierr);
  }

//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void _SUBR_(mvec_times_sdMat)(_ST_ alpha, _TYPE_(const_mvec_ptr) vV,
                                       _TYPE_(const_sdMat_ptr) vC,
                           _ST_ beta,  _TYPE_(mvec_ptr) vW,
                                       int* ierr)
  {
  *ierr=-99;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,C,vC,*ierr);
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
void _SUBR_(mvec_QR)(_TYPE_(const_mvec_ptr) V, 
        _TYPE_(mvec_ptr) Q, _TYPE_(sdMat_ptr) R, int* ierr)
  {
  // TODO - check the status of TSQR in Trilinos
  *ierr=-99;
  }

//!@}

} // extern "C"

