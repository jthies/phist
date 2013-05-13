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
void DcrsMat_read_mm(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  Epetra_CrsMatrix* A=NULL;
  *ierr=EpetraExt::MatrixMarketFileToCrsMatrix(filename,comm,A);
  *vA = (_TYPE_(crsMat_ptr))(A);
  }

//! read a matrix from a Ghost CRS (binary) file.
void DcrsMat_read_bin(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void DcrsMat_read_hb(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  *ierr=-99; // not implemented in epetra
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void DcrsMat_getrow_map(_TYPE_(const_crsMat_ptr) vA, map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const map_ptr_t)(&(A->RowMap()));
  }

//! get column distribution of a matrix
void DcrsMat_getcol_map(_TYPE_(const_crsMat_ptr) vA, map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const map_ptr_t)(&(A->ColMap()));
  }

//! get the map for vectors x in y=A*x
void DcrsMat_getdomain_map(_TYPE_(const_crsMat_ptr) vA, map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const map_ptr_t)(&(A->DomainMap()));
  }

//! get the map for vectors y in y=A*x
void DcrsMat_getrange_map(_TYPE_(const_crsMat_ptr) vA, map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const map_ptr_t)(&(A->RangeMap()));
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void Dmvec_create(const map_ptr_t vmap, int nvec, _TYPE_(mvec_ptr)* vV, 
        int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_Map, map,vmap,*ierr);
  Epetra_MultiVector* result = new Epetra_MultiVector(*map,nvec);
  *vV=(_TYPE_(mvec_ptr))(&result);
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void sdMat_create(int nrows, int ncols, _TYPE_(sdMat_ptr)* vM, double** val, int* ierr)
  {
  *ierr=0;
  Epetra_SerialComm comm;
  Epetra_LocalMap localMap(nrows,0,comm);
  Epetra_MultiVector* mv = new Epetra_MultiVector(localMap,ncols);
  if (mv==NULL) *ierr=-1;
  int lda;
  *ierr = mv->ExtractView(val,&lda);
  *vM=(_TYPE_(sdMat_ptr))mv;
  }

//@}

void Dmvec_extract_view(Dmvec_ptr_t vV, double** val, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V, vV, *ierr);
  int lda;
  *ierr = V->ExtractView(val,&lda);
  }

void DsdMat_extract_view(DsdMat_ptr_t vM, double** val, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,M, vM, *ierr);
  int lda;
  *ierr = M->ExtractView(val,&lda);
  }

//! \name destructors

//@{

//!
void DcrsMat_delete(_TYPE_(crsMat_ptr) vA, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_CrsMatrix,A,vA,*ierr);
  delete [] A;
  }

//!
void Dmvec_delete(_TYPE_(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  delete [] V;
  }

//!
void DsdMat_delete(_TYPE_(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,M,vM,*ierr);
  delete [] M;
  }

//@}

//! \name Numerical functions
//!@{

//! y=alpha*A*x+beta*y.
void DcrsMat_X_mvec(double alpha, _TYPE_(const_crsMat_ptr) vA, _TYPE_(const_mvec_ptr) vx, 
double beta, _TYPE_(const_mvec_ptr) vy, int* ierr)
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
void Dmvec_dot_mvec(_TYPE_(const_mvec_ptr) vV, _TYPE_(const_mvec_ptr) vW, double* s, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,W,vW,*ierr);
  _CHECK_ZERO_(V->Dot(*W,s),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
void DmvecT_X_mvec(double alpha, _TYPE_(const_mvec_ptr) vV, _TYPE_(const_mvec_ptr) vW, double beta, _TYPE_(sdMat_ptr) vC, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,C,vC,*ierr);
  _CHECK_ZERO_(C->Multiply('T','N',alpha,*V,*W,beta),*ierr);
  }


//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void Dmvec_X_sdMat(double alpha, _TYPE_(const_mvec_ptr) vV,
                                       _TYPE_(const_sdMat_ptr) vC,
                           _ST_ beta,  _TYPE_(mvec_ptr) vW,
                                       int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,C,vC,*ierr);
  _CHECK_ZERO_(W->Multiply('N', 'N', alpha, *v, *C, beta),*ierr);
  }
  
//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
void Dmvec_QR(_TYPE_(const_mvec_ptr) V, _TYPE_(mvec_ptr) Q, _TYPE_(sdMat_ptr) R, int* ierr)
  {
  // TODO - check the status of TSQR in Trilinos
  *ierr=-99;
  }

//!@}

