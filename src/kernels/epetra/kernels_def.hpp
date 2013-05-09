#include "Epetra_config.h"

#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "EpetraExt_CrsMatrixIn.h"

#include "../cpp_macros.h"

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void Dread_crsMat_mm(void** vA, char* filename,int* ierr)
  {
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  Epetra_CrsMatrix* A=NULL;
  *ierr=EpetraExt::MatrixMarketFileToCrsMatrix(filename,comm,A);
  *vA = (void*)(A);
  }

//! read a matrix from a Ghost CRS (binary) file.
void Dread_crsMat_bin(void** A, char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void Dread_crsMat_hb(void** vA, char* filename,int* ierr)
  {
  *ierr=-99; // not implemented in epetra
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void Dget_row_map(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (void*)(&(A->RowMap()));
  }

//! get column distribution of a matrix
void Dget_col_map(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (void*)(&(A->ColMap()));
  }

//! get the map for vectors x in y=A*x
void Dget_domain_map(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (void*)(&(A->DomainMap()));
  }

//! get the map for vectors y in y=A*x
void Dget_range_map(void* vA, void** vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (void*)(&(A->RangeMap()));
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void Dcreate_mvec(void* vmap, int nvec, void** vV, double** val, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_Map, map,vmap,*ierr);
  Epetra_MultiVector* result = new Epetra_MultiVector(*map,nvec);
  *vV=(void*)(&result);
  int lda;
  *ierr = result->ExtractView(val,&lda);
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void Dcreate_sdMat(int nrows, int ncols, void** vM, double** val, int* ierr)
  {
  *ierr=0;
  Epetra_SerialComm comm;
  Epetra_LocalMap localMap(nrows,0,comm);
  Epetra_MultiVector* mv = new Epetra_MultiVector(localMap,ncols);
  if (mv==NULL) *ierr=-1;
  int lda;
  *ierr = mv->ExtractView(val,&lda);
  *vM=(void*)mv;
  }

//@}

//! \name destructors

//@{

//!
void Ddelete_crsMat(void* vA, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_CrsMatrix,A,vA,*ierr);
  delete [] A;
  }

//!
void Ddelete_mvec(void* vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  delete [] V;
  }

//!
void Ddelete_sdMat(void* vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,M,vM,*ierr);
  delete [] M;
  }

//@}

//! \name Numerical functions
//!@{

//! y=alpha*A*x+beta*y.
void DcrsMat_X_mvec(double alpha, void* vA, void* vx, double beta, void* vy, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,x,vx,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,y,vy,*ierr);
  if (alpha==0.0)
    {
    if (beta==0.0)
      {
      _CHECK_ZERO_(y->PutScalar(0.0),*ierr);
      }
    else if (beta!=1.0)
      {
      _CHECK_ZERO_(y->Scale(beta),*ierr);
      }
    }
  else if (beta==0.0)
    {
    _CHECK_ZERO_(A->Multiply(false,*x,*y),*ierr);
    if (alpha!=1.0)
      {
      _CHECK_ZERO_(y->Scale(alpha),*ierr);
      }
    }
  else
    {
    Epetra_MultiVector Ax(y->Map(),y->NumVectors());
    _CHECK_ZERO_(A->Multiply(false,*x,Ax),*ierr);
    _CHECK_ZERO_(y->Update(alpha,Ax,beta),*ierr);
    }
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void Dmvec_dot_mvec(void* vv, void* vw, double* s, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,v,vv,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,w,vw,*ierr);
  _CHECK_ZERO_(v->Dot(*w,s),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
void DmvecT_X_mvec(double alpha, void* vV, void* vW, double beta, void* vC, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,C,vC,*ierr);
  _CHECK_ZERO_(C->Multiply('T','N',alpha,*V,*W,beta),*ierr);
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
void Dmvec_QR(void* V, void* Q, void* R, int* ierr)
  {
  // TODO - check the status of TSQR in Trilinos
  *ierr=-99;
  }

//!@}

