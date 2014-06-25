extern "C" {

// \name Matrix input from a file

//@{

// we implement only the double precision real type D
void SUBR(type_avail)(int* ierr)
  {
  *ierr=0;
  }

//! read a matrix from a MatrixMarket (ASCII) file
void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  Epetra_CrsMatrix* A=NULL;
  *ierr=EpetraExt::MatrixMarketFileToCrsMatrix(filename,comm,A);
  *vA = (TYPE(crsMat_ptr))(A);
  
/*  std::cerr << "filename was '"<<filename<<"'"<<std::endl;
  if (A==NULL) {std::cerr << "EpetraExt returned NULL"<<std::endl;}
  if (*ierr!=0) {std::cerr << "EpetraExt returned int code "<<*ierr<<std::endl;}
  */
  }

//! read a matrix from a Ghost CRS (binary) file.
void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  *ierr=-99; // not implemented in epetra
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->RowMap()));
  }

//! get column distribution of a matrix
void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->ColMap()));
  }

//! get the map for vectors x in y=A*x
void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->DomainMap()));
  }

//! get the map for vectors y in y=A*x
void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*ierr);
  *vmap = (const_map_ptr_t)(&(A->RangeMap()));
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
const_map_ptr_t vmap, int nvec, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_BlockMap, map,vmap,*ierr);
  Epetra_MultiVector* result;
  TRY_CATCH(result = new Epetra_MultiVector(*map,nvec),*ierr);
  if (result==NULL) *ierr=-1;
  *vV=(TYPE(mvec_ptr))(result);
  }

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values).
void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, const_map_ptr_t vmap, 
        _ST_* values, lidx_t lda, int nvec,
        int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_BlockMap, map,vmap,*ierr);
  Epetra_MultiVector* V = new Epetra_MultiVector(View, *map, values, lda, nvec);
  *vV=(TYPE(mvec_ptr))(V);
  return;
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        const_comm_ptr_t vcomm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_Comm, comm,vcomm,*ierr);
  Epetra_LocalMap localMap(nrows,0,*comm);
  Epetra_MultiVector* mv = new Epetra_MultiVector(localMap,ncols);
  if (mv==NULL) *ierr=-1;
  *vM=(TYPE(sdMat_ptr))mv;
  }

//@}

//! retrieve local length of the vectors in V
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) vV, lidx_t* len, int* ierr)
  {
  *ierr = 0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  *len = V->MyLength();
  }

//! retrieve the map of the vectors in V
void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  *vmap=(const_map_ptr_t)&V->Map();
  }

//! retrieve the comm used for MPI communication in V
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) vV, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  *vcomm=(const_comm_ptr_t)&V->Comm();
  }

//! retrieve number of vectors/columns in V
void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* ierr)
  {
  *ierr = 0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  *nvec = V->NumVectors();
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,M,vM,*ierr);
  *nrows = M->MyLength();
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,M,vM,*ierr);
  *ncols = M->NumVectors();
  }


void SUBR(mvec_extract_view)(Dmvec_ptr_t vV, double** val, lidx_t* lda, int* ierr)
  {
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V, vV, *ierr);
  PHIST_CHK_IERR(*ierr=V->ExtractView(val,lda),*ierr);
  }

void SUBR(sdMat_extract_view)(DsdMat_ptr_t vM, double** val, lidx_t* lda, int* ierr)
  {
  CAST_PTR_FROM_VOID(Epetra_MultiVector,M, vM, *ierr);
  PHIST_CHK_IERR(*ierr=M->ExtractView(val,lda),*ierr);
  }

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
                             TYPE(mvec_ptr)* vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  Epetra_MultiVector *Vblock;
  TRY_CATCH(Vblock = new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1),*ierr);
  if (*vVblock!=NULL)
    {
    CAST_PTR_FROM_VOID(Epetra_MultiVector,tmp,*vVblock,*ierr);
    delete tmp;
    }
  *vVblock = (TYPE(mvec_ptr))Vblock;
  }

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
                             TYPE(mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Vblock,vVblock,*ierr);
  Teuchos::RCP<const Epetra_MultiVector> Vcols;
  TRY_CATCH(Vcols = Teuchos::rcp(new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1)),*ierr);
  *Vblock = *Vcols;
  }

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,Vblock,vVblock,*ierr);
  Teuchos::RCP<Epetra_MultiVector> Vcols;
  TRY_CATCH(Vcols = Teuchos::rcp(new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1)),*ierr);
  *Vcols = *Vblock;
  }

//! not yet implemented
void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) vM,
                             TYPE(sdMat_ptr)* vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,M,vM,*ierr);
  Epetra_MultiVector *Mblock=NULL;
  if (*vMblock!=NULL)
    {
    CAST_PTR_FROM_VOID(Epetra_MultiVector,tmp,*vMblock,*ierr);
    delete tmp;
    }
  double* val;
  lidx_t lda;
  PHIST_CHK_IERR(*ierr=M->ExtractView(&val,&lda),*ierr);
  PHIST_DEB("in view (%d:%d,%d:%d), lda=%d\n",imin,imax,jmin,jmax,lda);
  Epetra_LocalMap localMap(imax-imin+1,M->Map().IndexBase(),M->Map().Comm());
  Mblock = new Epetra_MultiVector(View, localMap, val+imin+jmin*lda, lda, jmax-jmin+1);
  *vMblock = (TYPE(sdMat_ptr))Mblock;
  }

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM, 
                             TYPE(sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  TYPE(sdMat_ptr) vMtmp=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)((TYPE(sdMat_ptr))vM,&vMtmp,imin,imax,jmin,jmax,ierr),*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Mblock,vMblock,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Mtmp,vMtmp,*ierr);
  PHIST_CHK_IERR(*Mblock=*Mtmp,*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vMtmp,ierr),*ierr);
  }

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
                             TYPE(const_sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  TYPE(sdMat_ptr) vMtmp=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,&vMtmp,imin,imax,jmin,jmax,ierr),*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Mblock,vMblock,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Mtmp,vMtmp,*ierr);
  PHIST_CHK_IERR(*Mtmp=*Mblock,*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vMtmp,ierr),*ierr);
  }

//! \name destructors

//@{

//!
void SUBR(crsMat_delete)(TYPE(crsMat_ptr) vA, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_CrsMatrix,A,vA,*ierr);
  delete A;
  }

//!
void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  if (vV==NULL) return;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  delete V;
  }

//!
void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  if(vM==NULL) return;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,M,vM,*ierr);
  delete M;
  }

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, double value, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  PHIST_CHK_IERR(*ierr=V->PutScalar(value),*ierr);
  }

//! put scalar value into all elements of a multi-vector
void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vV, double value, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  PHIST_CHK_IERR(*ierr=V->PutScalar(value),*ierr);
  }

//! put random numbers into all elements of a multi-vector
void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  PHIST_CHK_IERR(*ierr=V->Random(),*ierr);
  }

void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* ierr)
  {
  *ierr = 0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  std::cout << *V;
  }

void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,M,vM,*ierr);
  std::cout << *M;
  }

//! put random numbers into all elements of a serial dense matrix
void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,M,vM,*ierr);
  PHIST_CHK_IERR(*ierr=M->Random(),*ierr);
  }

//! \name Numerical functions

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
  void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _ST_* vnrm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);  
  PHIST_CHK_IERR(*ierr=V->Norm2(vnrm),*ierr);
  return;
  }

  //! normalize (in the 2-norm) each column of v and return ||v||_2
  //! for each vector i in vnrm[i] (must be pre-allocated by caller)
  void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
                            _ST_* vnrm, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);  
  PHIST_CHK_IERR(*ierr=V->Norm2(vnrm),*ierr);
  for (int i=0;i<V->NumVectors();i++)
    {
    PHIST_CHK_IERR(*ierr=(*V)(i)->Scale(1.0/(_ST_)vnrm[i]),*ierr);
    }
  return;
  }

//! scale each column i of v and by scalar[i]
void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);  
  PHIST_CHK_IERR(*ierr=V->Scale(scalar),*ierr);
  return;
  }

//! scale each column i of v and by scalar[i]
void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            const _ST_* scalar, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);  
  for (int i=0;i<V->NumVectors();i++)
    {
    PHIST_CHK_IERR(*ierr=(*V)(i)->Scale(scalar[i]),*ierr);
    }
  return;
  }

//! y=alpha*x+beta*y
void SUBR(mvec_add_mvec)(double alpha, TYPE(const_mvec_ptr) vX,
                            double beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Y,vY,*ierr);
  PHIST_CHK_IERR(*ierr=Y->Update(alpha,*X,beta),*ierr);
  }

//! y[i]=alpha[i]*x[i]+beta*y[i], i=1..nvec
void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY,
                            int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,Y,vY,*ierr);

  for (int i=0;i<X->NumVectors(); i++)
    {
    TRY_CATCH((*Y)(i)->Update(alpha[i],*(*X)(i), beta),*ierr);
    }
  }


//! B=alpha*A+beta*B
void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,A,vA,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,B,vB,*ierr);
  PHIST_CHK_IERR(*ierr=B->Update(alpha,*A,beta),*ierr);  
  }

//! y=alpha*A*x+beta*y.
void SUBR(crsMat_times_mvec)(double alpha, TYPE(const_crsMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
double beta, TYPE(mvec_ptr) vy, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,x,vx,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,y,vy,*ierr);
  if (alpha==0.0)
    {
    if (beta==0.0)
      {
      PHIST_CHK_IERR(*ierr=y->PutScalar(0.0),*ierr);
      }
    else if (beta!=1.0)
      {
      PHIST_CHK_IERR(*ierr=y->Scale(beta),*ierr);
      }
    }
  else if (beta==0.0)
    {
    PHIST_CHK_IERR(*ierr=A->Multiply(false,*x,*y),*ierr);
    if (alpha!=1.0)
      {
      PHIST_CHK_IERR(*ierr=y->Scale(alpha),*ierr);
      }
    }
  else
    {
    Epetra_MultiVector Ax(y->Map(),y->NumVectors());
    PHIST_CHK_IERR(*ierr=A->Multiply(false,*x,Ax),*ierr);
    PHIST_CHK_IERR(*ierr=y->Update(alpha,Ax,beta),*ierr);
    }
  /*
  std::cout << *A << std::endl;
  std::cout << *((*x)(0)) << std::endl;
  std::cout << *((*y)(0)) << std::endl;
  */
  }

//! y=alpha*A*x+beta*y.
void SUBR(crsMatT_times_mvec)(double alpha, TYPE(const_crsMat_ptr) vA, TYPE(const_mvec_ptr) 
vx, 
double beta, TYPE(mvec_ptr) vy, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,x,vx,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,y,vy,*ierr);
  if (alpha==0.0)
    {
    if (beta==0.0)
      {
      PHIST_CHK_IERR(*ierr=y->PutScalar(0.0),*ierr);
      }
    else if (beta!=1.0)
      {
      PHIST_CHK_IERR(*ierr=y->Scale(beta),*ierr);
      }
    }
  else if (beta==0.0)
    {
    PHIST_CHK_IERR(*ierr=A->Multiply(true,*x,*y),*ierr);
    if (alpha!=1.0)
      {
      PHIST_CHK_IERR(*ierr=y->Scale(alpha),*ierr);
      }
    }
  else
    {
    Epetra_MultiVector Ax(y->Map(),y->NumVectors());
    PHIST_CHK_IERR(*ierr=A->Multiply(true,*x,Ax),*ierr);
    PHIST_CHK_IERR(*ierr=y->Update(alpha,Ax,beta),*ierr);
    }
  }

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
void SUBR(crsMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;

  PHIST_CHK_IERR(SUBR(crsMat_times_mvec)(alpha, A, x, beta, y, ierr), *ierr);
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &nvec, ierr), *ierr);
  _ST_ alpha_shifts[nvec];
  for(int i = 0; i < nvec; i++)
    alpha_shifts[i] = alpha*shifts[i];
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha_shifts, x, st::one(), y, ierr), *ierr);
}
//! dot product of vectors v_i and w_i, i=1..numvecs
void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, double* s, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*ierr);
  PHIST_CHK_IERR(*ierr=V->Dot(*W,s),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
void SUBR(mvecT_times_mvec)(double alpha, TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, double beta, TYPE(sdMat_ptr) vC, int* ierr)
  {
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*ierr);
  PHIST_CHK_IERR(*ierr=C->Multiply('T','N',alpha,*V,*W,beta),*ierr);
  }


//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void SUBR(mvec_times_sdMat)(double alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           double beta,  TYPE(mvec_ptr) vW,
                                       int* ierr)
  {
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,W,vW,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*ierr);
  PHIST_CHK_IERR(*ierr=W->Multiply('N', 'N', alpha, *V, *C, beta),*ierr);
  }

//! n x m serial dense matrix times m x k serial dense matrix gives n x k serial dense matrix,
//! C=alpha*V*W + beta*C
void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                           TYPE(const_sdMat_ptr) vW,
                               _ST_ beta, TYPE(sdMat_ptr) vC,
                                       int* ierr)
  {
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*ierr);
  PHIST_CHK_IERR(*ierr=C->Multiply('N', 'N', alpha, *V, *W, beta),*ierr);
  }


//! n x m transposed serial dense matrix times m x k serial dense matrix gives m x k serial dense matrix,
//! C=alpha*V*W + beta*C
void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                           TYPE(const_sdMat_ptr) vW,
                               _ST_ beta, TYPE(sdMat_ptr) vC,
                                       int* ierr)
  {
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*ierr);
  PHIST_CHK_IERR(*ierr=C->Multiply('T', 'N', alpha, *V, *W, beta),*ierr);
  }


//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, ierr>0   
//! indicates the dimension of the null-space of V. The first m-ierr    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.  
void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* ierr)
  {
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);

  *ierr=0;
  CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector,R,vR,*ierr);

  int rank;
  MT rankTol=32*mt::eps();
  if (V->NumVectors()==1)
    {
    // we need a special treatment here because TSQR
    // uses a relative tolerance to determine rank deficiency,
    // so a single zero vector is not detected to be rank deficient.
    MT nrm;
    PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,ierr),*ierr);
    PHIST_DEB("single vector QR, R=%8.4e\n",nrm);
    ST* Rval=NULL;
    int ldR;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(R,&Rval,&ldR,ierr),*ierr);
    rank=1;
    if (nrm<rankTol)
      {
      PHIST_DEB("zero vector detected\n");
      // randomize the vector
      PHIST_CHK_IERR(SUBR(mvec_random)(vV,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,ierr),*ierr);
      rank=0;// dimension of null space
      }
    *Rval=(ST)nrm;
    *ierr=1-rank;
    return;
    }

  if (R->ConstantStride()==false)
    {
    *ierr = -1;
    return;
    }
  int stride = R->Stride();
  int nrows = R->MyLength();
  int ncols = R->NumVectors();
    
#ifdef TESTING
  PHIST_CHK_IERR(*ierr=nrows-ncols,*ierr);
  PHIST_CHK_IERR(*ierr=nrows-V->NumVectors(),*ierr);
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
  params->set("randomizeNullSpace",true);
  tsqr.setParameterList(params);

  TRY_CATCH(rank = tsqr.normalize(*V,R_view),*ierr);  
  *ierr = ncols-rank;// return positive number if rank not full.
  return;
  }


//!@}

}// extern "C"
