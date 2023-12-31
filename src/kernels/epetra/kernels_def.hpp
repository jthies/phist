/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
// \name Matrix input from a file

//@{

// we implement only the double precision real type D
extern "C" void SUBR(type_avail)(int* iflag)
{
  *iflag=0;
}

//! read a matrix from a MatrixMarket (ASCII) file
extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* vA, phist_const_comm_ptr vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  bool repart = *iflag&PHIST_SPARSEMAT_PERM_GLOBAL;
  *iflag=0;
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Comm,comm,vcomm,*iflag);
  Epetra_CrsMatrix* A=NULL;
#if defined(EPETRA_NO_64BIT_GLOBAL_INDICES)||defined(PHIST_FORCE_32BIT_GIDX)
  *iflag=EpetraExt::MatrixMarketFileToCrsMatrix(filename,*comm,A);
#else
  *iflag=EpetraExt::MatrixMarketFileToCrsMatrix64(filename,*comm,A);
#endif

if (repart)
{
#ifdef PHIST_HAVE_ISORROPIA
    Teuchos::RCP<Epetra_Map> newRowMap=Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> newMatrix=Teuchos::null;
    PHIST_CHK_IERR(phist::epetra_internal::repartition(Teuchos::rcp(A,false),newRowMap,newMatrix,true,iflag),*iflag);
    delete A;
    A=newMatrix.release().getRawPtr();
#else
    PHIST_SOUT(PHIST_WARNING,"%s, not repartitioning the matrix because Isorropia is not available\n",__FUNCTION__);
#endif
}


  *vA = (TYPE(sparseMat_ptr))(A);
  
/*  std::cerr << "filename was '"<<filename<<"'"<<std::endl;
  if (A==NULL) {std::cerr << "EpetraExt returned NULL"<<std::endl;}
  if (*iflag!=0) {std::cerr << "EpetraExt returned int code "<<*iflag<<std::endl;}
  */
}

//! read a matrix from a MatrixMarket (ASCII) file
extern "C" void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* vA, phist_const_context_ptr vctx,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(const phist::internal::default_context,ctx,vctx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Map,map,ctx->row_map,*iflag);
  Epetra_CrsMatrix* A=NULL;
#if defined(EPETRA_NO_64BIT_GLOBAL_INDICES)||defined(PHIST_FORCE_32BIT_GIDX)
  *iflag=EpetraExt::MatrixMarketFileToCrsMatrix(filename,*map,A);
#else
  *iflag=EpetraExt::MatrixMarketFileToCrsMatrix64(filename,*map,A);
#endif

  *vA = (TYPE(sparseMat_ptr))(A);
}

//! read a matrix from a Ghost CRS (binary) file.
extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* vA, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // TODO - not implemented (should read the binary file format defined by ghost)
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! read a matrix from a Ghost CRS (binary) file.
extern "C" void SUBR(sparseMat_read_bin_with_context)(TYPE(sparseMat_ptr)* vA, phist_const_context_ptr vctx,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // TODO - not implemented (should read the binary file format defined by ghost)
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! read a matrix from a Harwell-Boeing (HB) file
extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* vA, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED; // not implemented in epetra
}

//! read a matrix from a Harwell-Boeing (HB) file
extern "C" void SUBR(sparseMat_read_hb_with_context)(TYPE(sparseMat_ptr)* vA, phist_const_context_ptr vctx,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED; // not implemented in epetra
}
//!@}

extern "C" void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA, phist_const_context_ptr vctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag)
{
  int iflag_in=*iflag;
  bool repart = (*iflag)&PHIST_SPARSEMAT_PERM_GLOBAL;
  bool ownMaps= (*iflag)&PHIST_SPARSEMAT_OWN_MAPS;
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(phist::internal::default_context,ctx,vctx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_Map,map,ctx->row_map,*iflag);
  PHIST_SOUT(PHIST_DEBUG,"sparseMat_create_fromRowFuncAndContext with iflag=%d\n",iflag_in);
  phist_gidx cols[maxnne];
  double vals[maxnne];

  // note: if the context contains a col_map, we do not force it onto the matrix anyway.
  //       this is because we don't have a fused_spmv_pair kernel for Epetra, so restricting
  //       the sparsity pattern of subsequently created matrices doesn't help.
 
  Epetra_CrsMatrix* A=NULL;
  PHIST_TRY_CATCH(A   = new Epetra_CrsMatrix(Copy,*map,maxnne),*iflag);
  for (phist_lidx i=0; i<A->NumMyRows(); i++)
  {
#if defined(EPETRA_NO_64BIT_GLOBAL_INDICES)||defined(PHIST_FORCE_32BIT_GIDX)
    ghost_gidx row = (ghost_gidx)map->GID(i);
#else
    ghost_gidx row = (ghost_gidx)map->GID64(i);
#endif
    ghost_lidx row_nnz;
    
    PHIST_CHK_IERR(*iflag=rowFunPtr(row,&row_nnz,cols,vals,last_arg),*iflag);
    PHIST_TRY_CATCH(A->InsertGlobalValues(row,row_nnz,vals,cols),*iflag);
  }
  if (ctx->range_map==NULL || ctx->domain_map==NULL)
  {
    PHIST_TRY_CATCH(A->FillComplete(),*iflag);
  }
  else
  {
    const Epetra_Map* domain_map=(const Epetra_Map*)(ctx->domain_map);
    const Epetra_Map* range_map=(const Epetra_Map*)(ctx->range_map);
    PHIST_TRY_CATCH(A->FillComplete(*domain_map,*range_map),*iflag);
  }

  *vA = (TYPE(sparseMat_ptr))(A);  

  if (repart)
  {
#ifdef PHIST_HAVE_ISORROPIA
    Teuchos::RCP<Epetra_Map> newRowMap=Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> dummy=Teuchos::null;
    PHIST_CHK_IERR(phist::epetra_internal::repartition(Teuchos::rcp(A,false),newRowMap,dummy,false,iflag),*iflag);
    if (newRowMap==Teuchos::null)
    {
      *iflag=PHIST_BAD_CAST;
      return;
    }
    if (ownMaps)
    {
      delete (phist::internal::default_context*)ctx;
    }
    phist::internal::default_context *new_ctx=new 
        phist::internal::default_context(newRowMap.get(),NULL,NULL);
    // create the matrix again and use this new map instead
    delete A;
    *iflag=iflag_in & ~PHIST_SPARSEMAT_PERM_GLOBAL;
    PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncAndContext)(vA,new_ctx,maxnne,rowFunPtr,last_arg,iflag),*iflag);
    A=(Epetra_CrsMatrix*)(*vA);
    delete new_ctx;

    // TODO: print statistics on the partitioning? Isorropia has nice tools for this in examples/

#else
    PHIST_SOUT(PHIST_WARNING,"matrix repartitioning with Epetra requires Isorropia (and Zoltan)\n");
#endif
  }
  else if (ownMaps)
  {
    ::phist::internal::contextCollection[*vA]=ctx;
  }
  return;
}

extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *vA, phist_const_comm_ptr vcomm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  int iflag_in=*iflag;
  *iflag=0;

  phist_map_ptr row_and_range_map=NULL;
  PHIST_CHK_IERR(phist_map_create(&row_and_range_map,vcomm,nrows,iflag),*iflag);
  phist_map_ptr domain_map=row_and_range_map;
  if (nrows!=ncols)
  {
    domain_map=NULL;
    PHIST_CHK_IERR(phist_map_create(&domain_map,vcomm,ncols,iflag),*iflag);
  }
  phist::internal::default_context *ctx=NULL;
  if (domain_map!=row_and_range_map)
  {
    ctx=new phist::internal::default_context(row_and_range_map,NULL,row_and_range_map,domain_map);
  }
  else
  {
    ctx=new phist::internal::default_context(row_and_range_map,NULL,NULL);
  }
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncAndContext)(vA,ctx,
        maxnne,rowFunPtr,last_arg,iflag),*iflag);
  if (domain_map!=row_and_range_map)
  {
    PHIST_CHK_IERR(phist_map_delete(domain_map,iflag),*iflag);
  }
  PHIST_CHK_IERR(phist_map_delete(row_and_range_map,iflag),*iflag);
  delete ctx;
  return;
}

//! \name get information about the data distribution in a matrix (maps)

//!@{
extern "C" void SUBR(sparseMat_local_nnz)(TYPE(const_sparseMat_ptr) vA, int64_t* local_nnz, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  *local_nnz = static_cast<int64_t>(A->NumMyNonzeros());
}

extern "C" void SUBR(sparseMat_global_nnz)(TYPE(const_sparseMat_ptr) vA, int64_t* global_nnz, int* iflag)
{
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  *global_nnz = static_cast<int64_t>(A->NumGlobalNonzeros());
}

//! get the row distribution of the matrix
extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(&(A->RowMap()));
}

//! get column distribution of a matrix
extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(&(A->ColMap()));
}

//! get the map for vectors x in y=A*x
extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(&(A->DomainMap()));
}

//! get the map for vectors y in y=A*x
extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(&(A->RangeMap()));
}
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
phist_const_map_ptr vmap, int nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_BlockMap, map,vmap,*iflag);
  Epetra_MultiVector* result;
  PHIST_TRY_CATCH(result = new Epetra_MultiVector(*map,nvec),*iflag);
  if (result==NULL) *iflag=-1;
  *vV=(TYPE(mvec_ptr))(result);
}

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        phist_const_comm_ptr vcomm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  Epetra_SerialComm serialComm;
  const Epetra_Comm* comm = NULL;
  if (vcomm!=NULL)
    comm=(const Epetra_Comm*)vcomm;
  else
    comm = &serialComm;
  Epetra_LocalMap localMap(nrows,0,*comm);
  Epetra_MultiVector* mv = new Epetra_MultiVector(localMap,ncols);
  if (mv==NULL) *iflag=-1;
  *vM=(TYPE(sdMat_ptr))mv;
}

//@}

//! retrieve the map of the vectors in V
extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  *vmap=(phist_const_map_ptr)&V->Map();
}

//! retrieve number of vectors/columns in V
extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  *nvec = V->NumVectors();
}

//! get number of cols in local dense matrix
extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,M,vM,*iflag);
  *nrows = M->MyLength();
}

//! get number of cols in local dense matrix
extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,M,vM,*iflag);
  *ncols = M->NumVectors();
}


extern "C" void SUBR(mvec_extract_view)(phist_Dmvec_ptr vV, double** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V, vV, *iflag);
  PHIST_CHK_IERR(*iflag=V->ExtractView(val,lda),*iflag);
}

extern "C" void SUBR(sdMat_extract_view)(phist_DsdMat_ptr vM, double** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,M, vM, *iflag);
  PHIST_CHK_IERR(*iflag=M->ExtractView(val,lda),*iflag);
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,in, v_in, *iflag);
  PHIST_CAST_PTR_FROM_VOID(      Epetra_MultiVector,out, v_out, *iflag);

  Teuchos::RCP<Epetra_Import> import = Teuchos::rcp(new Epetra_Import(out->Map(),in->Map()));
  *iflag=out->Import(*in,*import,Insert);
  return;
}

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
                             TYPE(mvec_ptr)* vVblock,
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  Epetra_MultiVector *Vblock;
  PHIST_TRY_CATCH(Vblock = new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1),*iflag);
  if (*vVblock!=NULL)
    {
    PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,tmp,*vVblock,*iflag);
    delete tmp;
    }
  *vVblock = (TYPE(mvec_ptr))Vblock;
}

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
                             TYPE(mvec_ptr) vVblock,
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Vblock,vVblock,*iflag);
  Teuchos::RCP<const Epetra_MultiVector> Vcols;
  PHIST_TRY_CATCH(Vcols = Teuchos::rcp(new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1)),*iflag);
  *Vblock = *Vcols;
}

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,Vblock,vVblock,*iflag);
  Teuchos::RCP<Epetra_MultiVector> Vcols;
  PHIST_TRY_CATCH(Vcols = Teuchos::rcp(new Epetra_MultiVector(View, *V, jmin, jmax-jmin+1)),*iflag);
  *Vcols = *Vblock;
}

//! not yet implemented
extern "C" void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) vM,
                             TYPE(sdMat_ptr)* vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,M,vM,*iflag);
  Epetra_MultiVector *Mblock=NULL;
  if (*vMblock!=NULL)
    {
    PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,tmp,*vMblock,*iflag);
    delete tmp;
    }
  double* val;
  phist_lidx lda;
  PHIST_CHK_IERR(*iflag=M->ExtractView(&val,&lda),*iflag);
  //PHIST_DEB("in view (%d:%d,%d:%d), lda=%d\n",imin,imax,jmin,jmax,lda);
  Epetra_LocalMap localMap(imax-imin+1,M->Map().IndexBase(),M->Map().Comm());
  Mblock = new Epetra_MultiVector(View, localMap, val+imin+jmin*lda, lda, jmax-jmin+1);
  *vMblock = (TYPE(sdMat_ptr))Mblock;
}

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM, 
                             TYPE(sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  TYPE(sdMat_ptr) vMtmp=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)((TYPE(sdMat_ptr))vM,&vMtmp,imin,imax,jmin,jmax,iflag),*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Mblock,vMblock,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Mtmp,vMtmp,*iflag);
  PHIST_CHK_IERR(*Mblock=*Mtmp,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vMtmp,iflag),*iflag);
}

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
                             TYPE(const_sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  TYPE(sdMat_ptr) vMtmp=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,&vMtmp,imin,imax,jmin,jmax,iflag),*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Mblock,vMblock,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Mtmp,vMtmp,*iflag);
  PHIST_CHK_IERR(*Mtmp=*Mblock,*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vMtmp,iflag),*iflag);
}

//! \name destructors

//@{

//!
extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) vA, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if(vA==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(Epetra_CrsMatrix,A,vA,*iflag);
  // this is to avoid memory leaks, the function sparseMat_get_context will create
  // a small wrapper object and store it in a map, associated with this pointer to
  // a sparseMat.
  phist::internal::delete_default_context(vA);
  delete A;
}

//!
extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (vV==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  delete V;
}

//!
extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if(vM==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,M,vM,*iflag);
  delete M;
}

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, double value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=V->PutScalar(value),*iflag);
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) vV,
        phist_mvec_elemFunc funPtr, void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  for (phist_lidx i=0;i<V->MyLength();i++)
  {
    for (int j=0; j<V->NumVectors(); j++)
    {
#if defined(EPETRA_NO_64BIT_GLOBAL_INDICES)||defined(PHIST_FORCE_32BIT_GIDX)
    phist_gidx row=V->Map().GID(i);
#else
    phist_gidx row=V->Map().GID64(i);
#endif
      PHIST_CHK_IERR(*iflag=funPtr(row,j,V->Pointers()[j]+i,last_arg),*iflag);
    }
  }
}

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vV, double value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=V->PutScalar(value),*iflag);
}

#ifndef PHIST_BUILTIN_RNG
//! put random numbers into all elements of a multi-vector
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=V->Random(),*iflag);
}
#endif

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  std::cout << *V;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,M,vM,*iflag);
  // may hang if called by just one MPI process
  //std::cout << *M;
  
  std::stringstream sos;
  for (int i=0; i< M->MyLength(); i++)
  {
    for (int j=0;j<M->NumVectors();j++)
    {
      sos << std::scientific 
          << std::setprecision(6)
          << std::setw(16) 
          << (*M)[j][i];
    }
    sos<<std::endl;
  }
  std::cout << sos.str();
}

#ifndef PHIST_BUILTIN_RNG
//! put random numbers into all elements of a serial dense matrix
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,M,vM,*iflag);
  PHIST_CHK_IERR(*iflag=M->Random(),*iflag);
}
#endif

//! put identity matrix into a small dense matrix \ingroup sdmat
extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) V, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;

  _ST_ *V_raw = NULL;
  phist_lidx lda;
  int m, n;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(V, &V_raw, &lda, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(V, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(V, &n, iflag), *iflag);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++)
      V_raw[lda*i+j] = (i==j) ? st::one() : st::zero();
}

//! \name Numerical functions

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
  extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _ST_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);  
  PHIST_CHK_IERR(*iflag=V->Norm2(vnrm),*iflag);
  return;
}

  //! normalize (in the 2-norm) each column of v and return ||v||_2
  //! for each vector i in vnrm[i] (must be pre-allocated by caller)
  extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
                            _ST_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);  
  PHIST_CHK_IERR(*iflag=V->Norm2(vnrm),*iflag);
  for (int i=0;i<V->NumVectors();i++)
    {
    PHIST_CHK_IERR(*iflag=(*V)(i)->Scale(1.0/(_ST_)vnrm[i]),*iflag);
    }
  return;
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);  
  PHIST_CHK_IERR(*iflag=V->Scale(scalar),*iflag);
  return;
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);  
  for (int i=0;i<V->NumVectors();i++)
    {
    PHIST_CHK_IERR(*iflag=(*V)(i)->Scale(scalar[i]),*iflag);
    }
  return;
}

//! y=alpha*x+beta*y
extern "C" void SUBR(mvec_add_mvec)(double alpha, TYPE(const_mvec_ptr) vX,
                            double beta,  TYPE(mvec_ptr)       vY, 
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Y,vY,*iflag);
  PHIST_CHK_IERR(*iflag=Y->Update(alpha,*X,beta),*iflag);
}

//! y[i]=alpha[i]*x[i]+beta*y[i], i=1..nvec
extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Y,vY,*iflag);

  for (int i=0;i<X->NumVectors(); i++)
    {
    PHIST_TRY_CATCH((*Y)(i)->Update(alpha[i],*(*X)(i), beta),*iflag);
    }
}

extern "C" void SUBR(mvec_times_mvec_elemwise)(_ST_ alpha, TYPE(const_mvec_ptr) vX,
                                                  TYPE(mvec_ptr)       vY, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,Y,vY,*iflag);
  PHIST_CHK_IERR(*iflag=Y->Multiply(alpha,*X,*Y,0.0),*iflag);  
}



//! B=alpha*A+beta*B
extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,B,vB,*iflag);
  PHIST_CHK_IERR(*iflag=B->Update(alpha,*A,beta),*iflag);  
}

//! B=alpha*A+beta*B
extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  // simple workaround
  TYPE(sdMat_ptr) I = NULL;
  int m = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vA,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&I,m,m,NULL,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_identity)(I,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(alpha,vA,I,beta,vB,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(I,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMat_times_mvec)(double alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
double beta, TYPE(mvec_ptr) vy, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;

  PHIST_COUNT_MATVECS(vx);

  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,y,vy,*iflag);
  if (alpha==0.0)
    {
    if (beta==0.0)
      {
      PHIST_CHK_IERR(*iflag=y->PutScalar(0.0),*iflag);
      }
    else if (beta!=1.0)
      {
      PHIST_CHK_IERR(*iflag=y->Scale(beta),*iflag);
      }
    }
  else if (beta==0.0)
    {
    PHIST_CHK_IERR(*iflag=A->Multiply(false,*x,*y),*iflag);
    if (alpha!=1.0)
      {
      PHIST_CHK_IERR(*iflag=y->Scale(alpha),*iflag);
      }
    }
  else
    {
    Epetra_MultiVector Ax(y->Map(),y->NumVectors());
    PHIST_CHK_IERR(*iflag=A->Multiply(false,*x,Ax),*iflag);
    PHIST_CHK_IERR(*iflag=y->Update(alpha,Ax,beta),*iflag);
    }
  /*
  std::cout << *A << std::endl;
  std::cout << *((*x)(0)) << std::endl;
  std::cout << *((*y)(0)) << std::endl;
  */
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMatT_times_mvec)(double alpha, TYPE(const_sparseMat_ptr) vA, 
TYPE(const_mvec_ptr) vx, 
double beta, TYPE(mvec_ptr) vy, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,y,vy,*iflag);
  if (alpha==0.0)
    {
    if (beta==0.0)
      {
      PHIST_CHK_IERR(*iflag=y->PutScalar(0.0),*iflag);
      }
    else if (beta!=1.0)
      {
      PHIST_CHK_IERR(*iflag=y->Scale(beta),*iflag);
      }
    }
  else if (beta==0.0)
    {
    PHIST_CHK_IERR(*iflag=A->Multiply(true,*x,*y),*iflag);
    if (alpha!=1.0)
      {
      PHIST_CHK_IERR(*iflag=y->Scale(alpha),*iflag);
      }
    }
  else
    {
    Epetra_MultiVector Ax(y->Map(),y->NumVectors());
    PHIST_CHK_IERR(*iflag=A->Multiply(true,*x,Ax),*iflag);
    PHIST_CHK_IERR(*iflag=y->Update(alpha,Ax,beta),*iflag);
    }
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;

  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &nvec, iflag), *iflag);

  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(alpha, A, x, beta, y, iflag), *iflag);
  _ST_ alpha_shifts[nvec];
  for(int i = 0; i < nvec; i++)
    alpha_shifts[i] = alpha*shifts[i];
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alpha_shifts, x, st::one(), y, iflag), *iflag);
}
//! dot product of vectors v_i and w_i, i=1..numvecs
extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, double* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*iflag);
  PHIST_CHK_IERR(*iflag=V->Dot(*W,s),*iflag);
}

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
extern "C" void SUBR(mvecT_times_mvec)(double alpha, TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, double beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*iflag);
  PHIST_CHK_IERR(*iflag=C->Multiply('T','N',alpha,*V,*W,beta),*iflag);
}


//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
extern "C" void SUBR(mvec_times_sdMat)(double alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           double beta,  TYPE(mvec_ptr) vW,
                                       int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*iflag);
  PHIST_CHK_IERR(*iflag=W->Multiply('N', 'N', alpha, *V, *C, beta),*iflag);
}

//! n x m serial dense matrix times m x k serial dense matrix gives n x k serial dense matrix,
//! C=alpha*V*W + beta*C
extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                           TYPE(const_sdMat_ptr) vW,
                               _ST_ beta, TYPE(sdMat_ptr) vC,
                                       int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*iflag);
  PHIST_CHK_IERR(*iflag=C->Multiply('N', 'N', alpha, *V, *W, beta),*iflag);
}


//! n x m transposed serial dense matrix times m x k serial dense matrix gives m x k serial dense matrix,
//! C=alpha*V*W + beta*C
extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                           TYPE(const_sdMat_ptr) vW,
                               _ST_ beta, TYPE(sdMat_ptr) vC,
                                       int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*iflag);
  PHIST_CHK_IERR(*iflag=C->Multiply('T', 'N', alpha, *V, *W, beta),*iflag);
}

//! n x m transposed serial dense matrix times m x k serial dense matrix gives m x k serial dense matrix,
//! C=alpha*V*W + beta*C
extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                                         _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_MultiVector,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,C,vC,*iflag);
  PHIST_CHK_IERR(*iflag=C->Multiply('N', 'T', alpha, *V, *W, beta),*iflag);
}



//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, iflag>0   
//! indicates the dimension of the null-space of V. The first m-iflag    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.  
extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector,R,vR,*iflag);

  int rank;
  MT rankTol=1000*mt::eps();
  if (V->NumVectors()==1)
    {
    // we need a special treatment here because TSQR
    // uses a relative tolerance to determine rank deficiency,
    // so a single zero vector is not detected to be rank deficient.
    MT nrm;
    PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,iflag),*iflag);
    PHIST_DEB("single vector QR, R=%8.4e\n",nrm);
    ST* Rval=NULL;
    int ldR;
    PHIST_CHK_IERR(SUBR(sdMat_extract_view)(vR,&Rval,&ldR,iflag),*iflag);
    rank=1;
    *Rval=(ST)nrm;
    if (nrm<rankTol)
      {
      PHIST_DEB("zero vector detected\n");
      // randomize the vector
      PHIST_CHK_IERR(SUBR(mvec_random)(vV,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,iflag),*iflag);
      rank=0;// dimension of null space
      *Rval=st::zero();
      }
    *iflag=1-rank;
    return;
    }

  if (R->ConstantStride()==false)
    {
    *iflag = -1;
    return;
    }
  int stride = R->Stride();
  int nrows = R->MyLength();
  int ncols = R->NumVectors();
    
#ifdef PHIST_TESTING
  PHIST_CHK_IERR(*iflag=nrows-ncols,*iflag);
  PHIST_CHK_IERR(*iflag=nrows-V->NumVectors(),*iflag);
#endif  

  Teuchos::RCP<Teuchos_sdMat_t> R_view
        = CreateTeuchosViewNonConst(Teuchos::rcp(R,false),iflag);
  if (*iflag) return;
#ifdef HAVE_BELOS_TSQR
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

  PHIST_TRY_CATCH(rank = tsqr.normalize(*V,R_view),*iflag);  
  *iflag = ncols-rank;// return positive number if rank not full.
#else
  *iflag=-99;
#endif
  return;
}


//!@}

