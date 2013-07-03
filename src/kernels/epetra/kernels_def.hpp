#include "phist_macros.h"
#include "phist_typedefs.h"
#include "phist_kernels.h"

#include "Epetra_config.h"

#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "EpetraExt_CrsMatrixIn.h"

#include "epetra_helpers.h"

#include "BelosEpetraAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"

// \name Matrix input from a file

//@{

extern "C" {

// we implement only the double precision real type D
void phist_Dtype_avail(int* ierr)
  {
  *ierr=0;
  }

//! read a matrix from a MatrixMarket (ASCII) file
void phist_DcrsMat_read_mm(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  Epetra_CrsMatrix* A=NULL;
  *ierr=EpetraExt::MatrixMarketFileToCrsMatrix(filename,comm,A);
  *vA = (_TYPE_(crsMat_ptr))(A);
  
/*  std::cerr << "filename was '"<<filename<<"'"<<std::endl;
  if (A==NULL) {std::cerr << "EpetraExt returned NULL"<<std::endl;}
  if (*ierr!=0) {std::cerr << "EpetraExt returned int code "<<*ierr<<std::endl;}
  */
  }

//! read a matrix from a Ghost CRS (binary) file.
void phist_DcrsMat_read_bin(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void phist_DcrsMat_read_hb(_TYPE_(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  *ierr=-99; // not implemented in epetra
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void phist_DcrsMat_get_row_map(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->RowMap()));
  }

//! get column distribution of a matrix
void phist_DcrsMat_get_col_map(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->ColMap()));
  }

//! get the map for vectors x in y=A*x
void phist_DcrsMat_get_domain_map(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->DomainMap()));
  }

//! get the map for vectors y in y=A*x
void phist_DcrsMat_get_range_map(_TYPE_(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->RangeMap()));
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void phist_Dmvec_create(_TYPE_(mvec_ptr)* vV, 
const_map_ptr_t vmap, int nvec, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_BlockMap, map,vmap,*ierr);
  Epetra_MultiVector* result;
  _TRY_CATCH_(result = new Epetra_MultiVector(*map,nvec),*ierr);
  if (result==NULL) *ierr=-1;
  *vV=(_TYPE_(mvec_ptr))(result);
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void phist_DsdMat_create(_TYPE_(sdMat_ptr)* vM, int nrows, int ncols, int* ierr)
  {
  *ierr=0;
  Epetra_SerialComm comm;
  Epetra_LocalMap localMap(nrows,0,comm);
  Epetra_MultiVector* mv = new Epetra_MultiVector(localMap,ncols);
  if (mv==NULL) *ierr=-1;
  *vM=(_TYPE_(sdMat_ptr))mv;
  }

//@}

//! retrieve local length of the vectors in V
void phist_Dmvec_my_length(_TYPE_(const_mvec_ptr) vV, lidx_t* len, int* ierr)
  {
  *ierr = 0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  *len = V->MyLength();
  }

//! retrieve the map of the vectors in V
void _SUBR_(mvec_get_map)(_TYPE_(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  *vmap=(const_map_ptr_t)&V->Map();
  }

//! retrieve number of vectors/columns in V
void phist_Dmvec_num_vectors(_TYPE_(const_mvec_ptr) vV, int* nvec, int* ierr)
  {
  *ierr = 0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  *nvec = V->NumVectors();
  }

void phist_Dmvec_extract_view(Dmvec_ptr_t vV, double** val, lidx_t* lda, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V, vV, *ierr);
  _CHECK_ZERO_((*V)(0)->ExtractView(val,lda),*ierr);
  }

void phist_DsdMat_extract_view(DsdMat_ptr_t vM, double** val, lidx_t* lda, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,M, vM, *ierr);
  _CHECK_ZERO_(M->ExtractView(val,lda),*ierr);
  }

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
void _SUBR_(mvec_view_block)(_TYPE_(mvec_ptr) vV,
                             _TYPE_(mvec_ptr)* vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  Epetra_MultiVector *Vblock;
  _TRY_CATCH_(Vblock = new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1),*ierr);
  *vVblock = (_TYPE_(mvec_ptr))Vblock;
  }

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void _SUBR_(mvec_get_block)(_TYPE_(const_mvec_ptr) vV,
                             _TYPE_(mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,Vblock,vVblock,*ierr);
  Teuchos::RCP<const Epetra_MultiVector> Vcols;
  _TRY_CATCH_(Vcols = Teuchos::rcp(new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1)),*ierr);
  *Vblock = *Vcols;
  }

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void _SUBR_(mvec_set_block)(_TYPE_(mvec_ptr) vV,
                             _TYPE_(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,Vblock,vVblock,*ierr);
  Teuchos::RCP<Epetra_MultiVector> Vcols;
  _TRY_CATCH_(Vcols = Teuchos::rcp(new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1)),*ierr);
  *Vcols = *Vblock;
  }

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void _SUBR_(sdMat_get_block)(_TYPE_(const_mvec_ptr) M, 
                             _TYPE_(mvec_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void _SUBR_(mvec_set_block)(_TYPE_(sdMat_ptr) M, 
                             _TYPE_(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! \name destructors

//@{

//!
void phist_DcrsMat_delete(_TYPE_(crsMat_ptr) vA, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_CrsMatrix,A,vA,*ierr);
  delete A;
  }

//!
void phist_Dmvec_delete(_TYPE_(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  delete V;
  }

//!
void phist_DsdMat_delete(_TYPE_(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,M,vM,*ierr);
  delete M;
  }

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
void _SUBR_(mvec_put_value)(_TYPE_(mvec_ptr) vV, double value, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  _CHECK_ZERO_(V->PutScalar(value),*ierr);
  }

//! put scalar value into all elements of a multi-vector
void _SUBR_(sdMat_put_value)(_TYPE_(sdMat_ptr) vV, double value, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  _CHECK_ZERO_(V->PutScalar(value),*ierr);
  }

//! put random numbers into all elements of a multi-vector
void _SUBR_(mvec_random)(_TYPE_(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  _CHECK_ZERO_(V->Random(),*ierr);
  }

//! put random numbers into all elements of a serial dense matrix
void _SUBR_(sdMat_random)(_TYPE_(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,M,vM,*ierr);
  _CHECK_ZERO_(M->Random(),*ierr);
  }

//! \name Numerical functions

  //! normalize (in the 2-norm) each column of v and return ||v||_2
  //! for each vector i in vnrm[i] (must be pre-allocated by caller)
  void _SUBR_(mvec_normalize)(_TYPE_(mvec_ptr) vV,
                            _ST_* vnrm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);  
  _CHECK_ZERO_(V->Norm2(vnrm),*ierr);
  for (int i=0;i<V->NumVectors();i++)
    {
    _CHECK_ZERO_(V(i)->Scale(_ONE_/(_ST_)vnrm[i]),*ierr);
    }
  return;
  }

//! scale each column i of v and by scalar[i]
void _SUBR_(mvec_scale)(_TYPE_(mvec_ptr) vV, 
                            _ST_* scalar, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);  
  for (int i=0;i<V->NumVectors();i++)
    {
    _CHECK_ZERO_(V(i)->Scale(scalar[i]),*ierr);
    }
  return;
  }

//! y=alpha*x+beta*y
void _SUBR_(mvec_add_mvec)(double alpha, _TYPE_(const_mvec_ptr) vX,
                            double beta,  _TYPE_(mvec_ptr)       vY, 
                            int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,X,vX,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,Y,vY,*ierr);
  _CHECK_ZERO_(Y->Update(alpha,*X,beta),*ierr);
  }


//! y=alpha*A*x+beta*y.
void phist_DcrsMat_times_mvec(double alpha, _TYPE_(const_crsMat_ptr) vA, _TYPE_(const_mvec_ptr) vx, 
double beta, _TYPE_(mvec_ptr) vy, int* ierr)
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
  /*
  std::cout << *A << std::endl;
  std::cout << *((*x)(0)) << std::endl;
  std::cout << *((*y)(0)) << std::endl;
  */
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void phist_Dmvec_dot_mvec(_TYPE_(const_mvec_ptr) vV, _TYPE_(const_mvec_ptr) vW, double* s, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,W,vW,*ierr);
  _CHECK_ZERO_(V->Dot(*W,s),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
void phist_DmvecT_times_mvec(double alpha, _TYPE_(const_mvec_ptr) vV, _TYPE_(const_mvec_ptr) vW, double beta, _TYPE_(sdMat_ptr) vC, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,C,vC,*ierr);
  _CHECK_ZERO_(C->Multiply('T','N',alpha,*V,*W,beta),*ierr);
  }


//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void phist_Dmvec_times_sdMat(double alpha, _TYPE_(const_mvec_ptr) vV,
                                       _TYPE_(const_sdMat_ptr) vC,
                           double beta,  _TYPE_(mvec_ptr) vW,
                                       int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,C,vC,*ierr);
  _CHECK_ZERO_(W->Multiply('N', 'N', alpha, *V, *C, beta),*ierr);
  }
  
//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
void phist_Dmvec_QR(_TYPE_(mvec_ptr) vV, _TYPE_(sdMat_ptr) vR, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Epetra_MultiVector,R,vR,*ierr);

  if (R->ConstantStride()==false)
    {
    *ierr = -1;
    return;
    }
  int stride = R->Stride();
  int nrows = R->MyLength();
  int ncols = R->NumVectors();
    
#ifdef TESTING
  _CHECK_ZERO_(nrows-ncols,*ierr);
  _CHECK_ZERO_(nrows-V->NumVectors(),*ierr);
#endif  

  Teuchos::RCP<Teuchos_sdMat_t> R_view
        = CreateTeuchosViewNonConst(Teuchos::rcp(R,false),ierr);
  if (*ierr) return;
  
  Belos::TsqrOrthoManager<double, Epetra_MultiVector> tsqr("phist/epetra");
  Teuchos::RCP<const Teuchos::ParameterList> valid_params = 
        tsqr.getValidParameters();
  // faster but numerically less robust settings:
  Teuchos::RCP<const Teuchos::ParameterList> fast_params = 
        tsqr.getFastParameters();
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp
        (new Teuchos::ParameterList(*valid_params));
  tsqr.setParameterList(params);

  int rank;
  _TRY_CATCH_(rank = tsqr.normalize(*V,R_view),*ierr);  
  *ierr = ncols-rank;// return positive number if rank not full.
  return;
  }


//!@}

}// extern "C"
