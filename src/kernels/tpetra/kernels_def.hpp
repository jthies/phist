/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "Kokkos_View.hpp"
#include "Tpetra_MultiVector.hpp"

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
#if defined(IS_DOUBLE) && !defined(IS_COMPLEX) && defined(HAVE_TPETRA_INST_DOUBLE)
  *iflag = PHIST_SUCCESS;
#elif defined(IS_DOUBLE) && defined(IS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
  *iflag = PHIST_SUCCESS;
#elif !defined(IS_DOUBLE) && !defined(IS_COMPLEX) && defined(HAVE_TPETRA_INST_FLOAT)
  *iflag = PHIST_SUCCESS;
#elif !defined(IS_DOUBLE) && defined(IS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
  *iflag = PHIST_SUCCESS;
#endif
  return;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* matrixPtr,
                                        phist_const_comm_ptr vcomm,
                                        const char* fileName, int* iflag)
{ 
    PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
    if (fileName == nullptr)
    {
        *iflag = PHIST_INVALID_INPUT;
        return;
    }
    Tpetra::MatrixMarket::Reader<Traits<_ST_>::sparseMat_t> reader{};

    Teuchos::RCP<Traits<_ST_>::sparseMat_t> resultMatrix{};
    PHIST_CAST_PTR_FROM_VOID(const comm_type, comm, vcomm, *iflag);

    auto phist_comm_ptr = Teuchos::rcp(comm, false);

    PHIST_TRY_CATCH
    (
        resultMatrix = reader.readSparseFile(std::string{fileName}, phist_comm_ptr),
        *iflag       
    );

    auto resultMatrixPtr = resultMatrix.release();

    *matrixPtr = (TYPE(sparseMat_ptr))(resultMatrixPtr.get()); 
    *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, 
                                         phist_const_comm_ptr comm,
                                         const char* filename, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, 
                                        phist_const_comm_ptr comm,
                                        const char* filename, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* A, 
                                                     phist_const_context_ptr ctx,
                                                     const char* filename, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_bin_with_context)(TYPE(sparseMat_ptr)* A, 
                                                      phist_const_context_ptr ctx,
                                                      const char* filename, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb_with_context)(TYPE(sparseMat_ptr)* A, 
                                                     phist_const_context_ptr ctx,
                                                     const char* filename, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}

// Fill the sparse matrix row by row according to the context
extern "C" void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA, 
                                                             phist_const_context_ptr vctx,
                                                             phist_lidx maxnne, 
                                                             phist_sparseMat_rowFunc rowFunPtr, 
                                                             void* last_arg, int *iflag)
try
{
  int iflag_in = *iflag;
  PHIST_CAST_PTR_FROM_VOID(const phist::internal::default_context, ctx, vctx, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const map_type, tpetraMap, ctx->row_map, *iflag);

  // Check if the sparse mat may take possesion of the map
  bool ownMap = iflag_in & PHIST_SPARSEMAT_OWN_MAPS;
  auto mapPtr = Teuchos::rcp(tpetraMap, ownMap);
  
  auto sparseMat = new Traits<_ST_>::sparseMat_t(mapPtr, static_cast<int>(maxnne));

  phist_lidx numRows = sparseMat->getNodeNumRows();

  // note: parallel_for here leads to segfault, and I could not find a Tpetra example
  // where they insertGlobalVAlues in a parallel_for. Probably the function is not thread-safe
  // unless HAVE_TEUCHOS_THREADSAFE is defined.
  for (phist_lidx idx=0; idx<numRows; idx++)

  {
    *iflag = PHIST_SUCCESS;
    int iflag_local = PHIST_SUCCESS;
    phist_gidx cols[maxnne];
    _ST_ vals[maxnne];
    ghost_gidx row = tpetraMap->getGlobalElement(idx);
    ghost_lidx row_nnz;

    iflag_local = rowFunPtr(row, &row_nnz, cols, vals, last_arg);
    if (iflag_local) 
    {
      throw iflag_local;
    }

    if (row_nnz != 0)
    {
      Teuchos::ArrayView<phist_gidx> cols_v{cols, row_nnz};
      Teuchos::ArrayView<_ST_> vals_v{vals, row_nnz};
      
      sparseMat->insertGlobalValues(row, cols_v, vals_v);
    } 
  }

  const auto range_map = (const phist::tpetra::map_type*)(ctx->range_map);
  const auto domain_map = (const phist::tpetra::map_type*)(ctx->domain_map);

  if (range_map != nullptr && domain_map != nullptr)
  {
    auto range = Teuchos::rcp(range_map, ownMap);
    auto domain = Teuchos::rcp(domain_map, ownMap);
    PHIST_TRY_CATCH(sparseMat->fillComplete(domain, range),*iflag);
  }
  else
  {
    PHIST_TRY_CATCH(sparseMat->fillComplete(),*iflag);
  }

  *vA = (TYPE(sparseMat_ptr))(sparseMat);

  if (ownMap)
  {
    phist::internal::contextCollection[*vA]=(phist::internal::default_context*)ctx;
  }
  
}
catch (std::exception &ex)
{
  #ifdef PHIST_TESTING
    std::cout << "Caught exception: " << ex.what() << '\n';
  #endif
  *iflag = PHIST_CAUGHT_EXCEPTION;
}  
catch (...) 
{
  *iflag = PHIST_CAUGHT_EXCEPTION; 
}

extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, 
                                                   phist_const_comm_ptr comm,
                                                   phist_gidx nrows, 
                                                   phist_gidx ncols, 
                                                   phist_lidx maxnne,
                                                   phist_sparseMat_rowFunc rowFunPtr, 
                                                   void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  int iflag_in = *iflag;

  phist_map_ptr vmap = nullptr;
  PHIST_CHK_IERR(phist_map_create(&vmap, comm, nrows, iflag), *iflag);
  // we have to pass in a context object, but only the row map is actually needed to create the matrix:
  auto ctx = new phist::internal::default_context(vmap, nullptr, nullptr);
  //The matrix will take ownership of the map and context:
  *iflag = iflag_in | PHIST_SPARSEMAT_OWN_MAPS;

  // We delegate the construction of the matrix
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncAndContext)(A, ctx, maxnne, 
                                                              rowFunPtr, last_arg, 
                                                              iflag),
                 *iflag);
}

extern "C" void SUBR(sparseMat_local_nnz)(TYPE(const_sparseMat_ptr) A,
                                            int64_t* local_nnz, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  // note: I don't know how to get this property from Tpetra, but since we want to allow perfcheck
  // anyway I return iflag=0.
  *iflag=0;
  *local_nnz=0;
}                                            

extern "C" void SUBR(sparseMat_global_nnz)(TYPE(const_sparseMat_ptr) vA, int64_t* global_nnz, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  // note: I don't know how to get this property from Tpetra, but since we want to allow perfcheck
  // anyway I return iflag=0.
  *iflag=0;
  *global_nnz=0;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, 
                                            phist_const_map_ptr* map, 
                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);

  *map = (phist_const_map_ptr)(mat->getRowMap().get());

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, 
                                            phist_const_map_ptr* map, 
                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getColMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, 
                                               phist_const_map_ptr* map, 
                                               int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getDomainMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, 
                                              phist_const_map_ptr* map, 
                                              int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getRangeMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vec, 
                                  phist_const_map_ptr map, 
                                  phist_lidx nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const map_type, inMap, map, *iflag);
  PHIST_PERFCHECK_VERIFY_MVEC_CREATE(map,nvec,iflag);  
  auto map_ptr = Teuchos::rcp(inMap, false);
  auto result = new Traits<_ST_>::mvec_t(map_ptr, nvec);
  *vec = (TYPE(mvec_ptr))(result);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* mat, int nrows, int ncols, 
                                   phist_const_comm_ptr comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  const comm_type* localComm = (const comm_type*)comm;

  auto phistCommPtr = localComm == nullptr ? 
        Teuchos::DefaultComm<int>::getDefaultSerialComm(Teuchos::null)
      :
        Teuchos::rcp(localComm, false);

  auto localMap = Teuchos::rcp(new sdMat_map_type(nrows, 0, phistCommPtr, Tpetra::LocallyReplicated));
  auto matrix = new Traits<_ST_>::sdMat_t(localMap, ncols);

  *mat = (TYPE(sdMat_ptr))(matrix);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, phist_const_map_ptr* map, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mVec, V, *iflag);
  *map = (phist_const_map_ptr)(mVec->getMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mVec, V, *iflag);
  *nvec = mVec->getNumVectors();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMat, M, *iflag);
  *nrows = sdMat->getLocalLength();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMat, M, *iflag);
  *ncols = sdMat->getNumVectors();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_sourceVec, 
                                   TYPE(mvec_ptr) v_targetVec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, sourceVec, v_sourceVec, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, targetVec, v_targetVec, *iflag);

  auto import = Teuchos::rcp(new import_type(targetVec->getMap(), sourceVec->getMap()));

  PHIST_TRY_CATCH(targetVec->doImport(*sourceVec, *import, Tpetra::INSERT), *iflag);

  *iflag = PHIST_SUCCESS;
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mVec, V, *iflag);

  Teuchos::ArrayRCP<_ST_> valptr = mVec->get1dViewNonConst();
  *val = valptr.get();
  *lda = mVec->getStride();
  PHIST_CHK_IERR(*iflag = (*val == nullptr) ? PHIST_BAD_CAST: PHIST_SUCCESS, *iflag);
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMat, V, *iflag);

  Teuchos::ArrayRCP<_ST_> valptr = sdMat->get1dViewNonConst();
  *val = valptr.release();
  *lda = sdMat->getStride();

  PHIST_CHK_IERR(*iflag= (val==nullptr) ? PHIST_BAD_CAST : PHIST_SUCCESS, *iflag);
}

#ifdef PHIST_HIGH_PRECISION_KERNELS
extern "C" void SUBR(sdMat_extract_error)(TYPE(sdMat_ptr) M, _ST_** MC_raw, int* iflag)
{
  *iflag = PHIST_NOT_IMPLEMENTED;
}
#endif

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
                                      TYPE(mvec_ptr)* Vblock,
                                      int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, V, *iflag);

  Teuchos::RCP<Traits<_ST_>::mvec_t > block;
  PHIST_TRY_CATCH(block = mvec->subViewNonConst(Teuchos::Range1D(jmin, jmax)), *iflag);

  // check if Vblock already has data and if so delete it
  if (*Vblock != nullptr)
  {
    PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, tmp, *Vblock, *iflag);
    delete tmp;
  }
  *Vblock = (TYPE(mvec_ptr))(block.release().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
                                     TYPE(mvec_ptr) Vblock,
                                     int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecIn, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvecOut, Vblock, *iflag);

  Teuchos::RCP<const Traits<_ST_>::mvec_t> cols;
  PHIST_TRY_CATCH(cols = mvecIn->subView(Teuchos::Range1D(jmin, jmax)), *iflag);

  PHIST_TRY_CATCH(Tpetra::deep_copy(*mvecOut, *cols), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
                                     TYPE(const_mvec_ptr) Vblock,
                                     int jmin, int jmax, int* iflag) 
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvecDest, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecSrc, Vblock, *iflag);

  Teuchos::RCP<Traits<_ST_>::mvec_t> cols;
  PHIST_TRY_CATCH(cols = mvecDest->subViewNonConst(Teuchos::Range1D(jmin, jmax)), *iflag);

  PHIST_TRY_CATCH(Tpetra::deep_copy(*cols, *mvecSrc), *iflag);

  *iflag = PHIST_SUCCESS;  
}

extern "C" void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) vM, 
                                       TYPE(sdMat_ptr)* vMblock,
                                       int imin, int imax, 
                                       int jmin, int jmax, int* iflag)
{
  using std::make_pair;
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, mat, vM, *iflag);
  PHIST_PERFCHECK_VERIFY_SMALL;
  if (*vMblock != nullptr)
  {
    PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,tmp,*vMblock,*iflag);
    delete tmp;
  }

  // Get a view of the whole matrix
  auto view = mat->getDualView();

  // Take a subview, subview = view[imin:imax ; jmin: jmax]
  // Phists wants inclusive endpoints, so we correct by adding + 1
  auto subview = Kokkos::subview(view, make_pair(imin, imax + 1), 
                                       make_pair(jmin, jmax + 1));

  // create new local map
  auto map = Teuchos::rcp(new sdMat_map_type(imax - imin + 1, 0, 
                                       mat->getMap()->getComm(),
                                       Tpetra::LocallyReplicated));

  // Create sdMat with this map, which views the requested rows and cols
  Traits<_ST_>::sdMat_t* block = new Traits<_ST_>::sdMat_t(map, subview,view);

  *vMblock = (TYPE(sdMat_ptr))block;
  *iflag = PHIST_SUCCESS;
}

// get a new matrix that is a copy of some rows and columns of the original one,  
// Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
// and the corresponding columns of M are copied into the value array    
// of Mblock. M is not modified.
extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM, 
                             TYPE(sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  *iflag=0;
  TYPE(sdMat_ptr) vMview=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)((TYPE(sdMat_ptr))vM,&vMview,imin,imax,jmin,jmax,iflag),*iflag);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,Mview,vMview,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  PHIST_TRY_CATCH(Tpetra::deep_copy(*Mblock,*Mview),*iflag); // copy operation
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vMview,iflag),*iflag);
}

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
                             TYPE(const_sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  *iflag = PHIST_SUCCESS;
  TYPE(sdMat_ptr) vMview = nullptr;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM, &vMview, imin, imax, jmin, jmax, iflag), *iflag);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, Mblock, vMblock, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, Mview, vMview, *iflag);

  PHIST_TRY_CATCH(Tpetra::deep_copy(*Mview, *Mblock), *iflag); // copy operation

  PHIST_CHK_IERR(SUBR(sdMat_delete)(vMview, iflag), *iflag);
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // this is to avoid memory leaks, the function sparseMat_get_context will create
  // a small wrapper object and store it in a map, associated with this pointer to
  // a sparseMat.
  phist::internal::delete_default_context(A);
  if (A == nullptr) 
  {
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  delete mat;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) mvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  if (mvec == nullptr)
  {
    *iflag = PHIST_SUCCESS;
    return;
  }

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mVec, mvec, *iflag);
  delete mVec;

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) sdmat, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  if (sdmat == nullptr)
  {
    *iflag = PHIST_SUCCESS;
    return;
  }

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMat, sdmat, *iflag);
  delete sdMat;

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vec, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(vec,iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  PHIST_TRY_CATCH(mvec->putScalar(value), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
                                    phist_mvec_elemFunc funPtr, 
                                    void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, V, *iflag);

  auto data = mvec->get2dViewNonConst();
  phist_lidx numRow = mvec->getLocalLength();
  phist_lidx numCol = mvec->getNumVectors();
  for (int col = 0; col != numCol; ++col)
  {
    for (int row = 0; row != numRow; ++row)
    {
      phist_gidx globalRowIdx = mvec->getMap()->getGlobalElement(row);
      PHIST_CHK_IERR(*iflag = funPtr(globalRowIdx, col, 
                                     (void*)(&(data[col][row])),
                                     last_arg), 
                     *iflag);
    }
  }
}


extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) sdMat, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdmat, sdMat, *iflag);

  PHIST_TRY_CATCH(sdmat->putScalar(value), *iflag);

  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, V, *iflag);

  PHIST_TRY_CATCH(V->randomize(), *iflag);

  *iflag = PHIST_SUCCESS;
}
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  #include "phist_std_typedefs.hpp"  
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, mat, M, *iflag);

  #ifdef PHIST_HAVE_MPI
    // generate the same data on all processes,
    // by using an allreduction^^ not very nice, but works
    // TODO: improve this
    int myRank = M->getMap()->getComm()->getRank();
    if( myRank == 0 )
    {
      PHIST_TRY_CATCH(M->randomize(), *iflag);
    }
    else
    {
      PHIST_TRY_CATCH(M->putScalar(st::zero()), *iflag);
    }
    PHIST_TRY_CATCH(M->reduce(), *iflag);
  #else
    PHIST_TRY_CATCH(M->randomize(), *iflag);
  #endif
  *iflag = PHIST_SUCCESS;
}
#endif

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vec,*iflag);
  Teuchos::FancyOStream fos{Teuchos::rcp(&std::cout,false)};
  fos << std::scientific << std::setw(16) << std::setprecision(12);
  V->describe(fos,Teuchos::VERB_EXTREME);
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  Teuchos::FancyOStream fos{Teuchos::rcp(&std::cout,false)};
  fos << std::scientific << std::setw(12) << std::setprecision(6);
  // this hangs if the function is called by not all MPI ranks (see #108)
  M->describe(fos,Teuchos::VERB_EXTREME); 
}

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) mat, int* iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;

  _ST_* raw_values = nullptr; 

  phist_lidx lda;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(mat, &raw_values, &lda, iflag), *iflag);

  int numRows;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &numRows, iflag), *iflag);

  int numCols;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &numCols, iflag), *iflag);

  // Put ones on the diagonal and zeroes on the off diagonals
  for (int col = 0; col != numCols; ++col)
  {
    for (int row = 0; row != numRows; ++row)
    { // Branch predictor should do well enough here, [[unlikely]] attribute in c++20?
      raw_values[lda * col + row] = col == row ? st::one() : st::zero();
    }
  }

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vec,
                                 _MT_* vnrm, int* iflag)
{ 
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(vec,vec,iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  int nvec = mvec->getNumVectors();
  Teuchos::ArrayView<_MT_> norms{vnrm, nvec};

  PHIST_TRY_CATCH(mvec->norm2(norms), *iflag);

  *iflag = PHIST_SUCCESS;
}
// TODO: All the functions that call Tpetra::MultiVector::scale()
//       should use Kokkos views instead of Teuchos ArrayViews
extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vec,
                                     _MT_* vnrm, int* iflag)
{ 
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  int nvec = mvec->getNumVectors();
  Teuchos::ArrayView<_MT_> norms{vnrm, nvec};
  Teuchos::Array<_ST_> scaling{nvec};

  PHIST_TRY_CATCH(mvec->norm2(norms), *iflag);

  for (int idx = 0; idx != nvec; ++idx)
    scaling[idx] = norms[idx] == mt::zero() ? mt::one() : mt::one() / norms[idx];
  
  PHIST_TRY_CATCH(mvec->scale(scaling), *iflag);

  *iflag = PHIST_SUCCESS; 
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vec, 
                                 _ST_ scalar, int* iflag)
{  
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_SCALE(vec,iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);
  PHIST_TRY_CATCH(mvec->scale(scalar), *iflag);
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vec, 
                                  const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_SCALE(vec,iflag);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  int nvec = mvec->getNumVectors();
  Teuchos::ArrayView<_ST_> scalars{(_ST_*)scalar, nvec};

  PHIST_TRY_CATCH(mvec->scale(scalars), *iflag);

  *iflag = PHIST_SUCCESS; 
}

// vecOut = alpha * vecIn + beta * vecOut
extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vecIn,
                                    _ST_ beta, TYPE(mvec_ptr) vecOut, 
                                    int* iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(alpha,vecIn,beta,vecOut,iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, vec, vecIn, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, resultVec, vecOut, *iflag);

  // Check if we only want to copy vec into resultVec
  if (alpha == st::one() and beta == st::zero())
  {
    PHIST_TRY_CATCH(Tpetra::deep_copy(*resultVec, *vec), *iflag);
  }
  else
    PHIST_TRY_CATCH(resultVec->update(alpha, *vec, beta), *iflag);

  *iflag = PHIST_SUCCESS;
}

// vecOut[i] = alpha[i] * vecIn[i] + beta * vecOut[i]
extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vecIn,
                                     _ST_ beta, TYPE(mvec_ptr) vecOut, 
                                     int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(*alpha,vecIn,beta,vecOut,iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, vec, vecIn, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, resultVec, vecOut, *iflag);

  for (unsigned int idx = 0; idx != vec->getNumVectors(); ++idx)
  {
    PHIST_TRY_CATCH(resultVec->getVectorNonConst(idx)->update(alpha[idx], 
                                                           *vec->getVector(idx), 
                                                           beta),
                    *iflag);
  }
  *iflag = PHIST_SUCCESS;
}

// 
extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) matIn,
                                      _ST_ beta,  TYPE(sdMat_ptr) matOut, 
                                      int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn, matIn, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMatOut, matOut, *iflag);

  PHIST_TRY_CATCH(sdMatOut->update(alpha, *sdMatIn, beta), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                                       _ST_ beta,  TYPE(sdMat_ptr) B, 
                                       int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;

  // This cannot be easily done in Tpetra, so we will use a workaround
  // alpha * X^T + beta * Y = alpha * X^T * I + beta * Y
  TYPE(sdMat_ptr) identity = nullptr;
  int numCols;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(A, &numCols, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&identity, numCols, numCols, nullptr, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_identity)(identity, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(alpha, A, identity, beta, B, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(identity, iflag), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) A, 
                                                       TYPE(const_mvec_ptr) x, 
                                                       int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
                                           TYPE(const_mvec_ptr) x, _ST_ beta, 
                                           TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_COUNT_MATVECS(x);
  PHIST_PERFCHECK_VERIFY_SPMV(alpha,A,beta,x,y,0.0,0.0,iflag);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, matrix, A, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecIn, x, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvecOut, y, *iflag);

  PHIST_TRY_CATCH(matrix->apply(*mvecIn, *mvecOut, Teuchos::NO_TRANS, alpha, beta),
                  *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
                                            TYPE(const_mvec_ptr) x, _ST_ beta, 
                                            TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_COUNT_MATVECS(x);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, matrix, A, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecIn, x, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvecOut, y, *iflag);

  #ifdef IS_COMPLEX
    PHIST_TRY_CATCH(matrix->apply(*mvecIn, *mvecOut, Teuchos::CONJ_TRANS, alpha, beta),
                    *iflag);
  #else
    PHIST_TRY_CATCH(matrix->apply(*mvecIn, *mvecOut, Teuchos::TRANS, alpha, beta),
                    *iflag);
  #endif

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
                                                     const _ST_ shifts[], TYPE(const_mvec_ptr) x, 
                                                     _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // First compute the mvec_t y = alpha * A * x + beta * y
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(alpha, A, x, beta, y, iflag),
                 *iflag);

  int numVec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &numVec, iflag),
                 *iflag);

  _ST_ alphaShifts[numVec];
  for(int idx = 0; idx != numVec; ++idx)
  {
    alphaShifts[idx] = alpha * shifts[idx];
  }
  // Add the shifts to y
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alphaShifts, x, st::one(), y, iflag),
                 *iflag);

  *iflag = PHIST_SUCCESS;
}


extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
                                    TYPE(const_mvec_ptr) w, 
                                    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(v,w,iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvec1, v, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvec2, w, *iflag);

  Teuchos::ArrayView<_ST_> result{s, (int)mvec1->getNumVectors()};

  PHIST_TRY_CATCH(mvec1->dot(*mvec2, result), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_sdMat_ptr) C, 
                                       _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(alpha,V,beta,W,iflag);
  
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecIn, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdmat, C, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvecOut, W, *iflag);

  PHIST_TRY_CATCH(mvecOut->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                                    alpha, *mvecIn, *sdmat, beta),
                  *iflag);

  //PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  *iflag = PHIST_SUCCESS;
}


extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                        TYPE(const_sdMat_ptr) W, 
                                        _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn1, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn2, W, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMatOut, C, *iflag);

  PHIST_TRY_CATCH(sdMatOut->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                                    alpha, *sdMatIn1, *sdMatIn2, beta),
                  *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                         TYPE(const_sdMat_ptr) W, 
                                         _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn1, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn2, W, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMatOut, C, *iflag);

  PHIST_TRY_CATCH(sdMatOut->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
                                    alpha, *sdMatIn1, *sdMatIn2, beta),
                  *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                         TYPE(const_sdMat_ptr) W, 
                                         _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn1, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMatIn2, W, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMatOut, C, *iflag);

  PHIST_TRY_CATCH(sdMatOut->multiply(Teuchos::NO_TRANS, Teuchos::CONJ_TRANS,
                                    alpha, *sdMatIn1, *sdMatIn2, beta),
                  *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
                                       TYPE(const_mvec_ptr) W, 
                                       _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag);
  
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecIn1, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvecIn2, W, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvecOut, C, *iflag);

  PHIST_TRY_CATCH(mvecOut->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
                                    alpha, *mvecIn1, *mvecIn2, beta),
                  *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Dmvec* reV, phist_Dmvec* imV, int *iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvec, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_MT_>::mvec_t, reMvec, reV, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_MT_>::mvec_t, imMvec, imV, *iflag);
  
  //mvec->sync<Kokkos::HostSpace> ();
  reMvec->sync<Kokkos::HostSpace> ();
  imMvec->sync<Kokkos::HostSpace> ();

  auto const sourceMvecView = mvec->getLocalView<Kokkos::HostSpace>();
  auto reMvecView = reMvec->getLocalView<Kokkos::HostSpace>();
  auto imMvecView = imMvec->getLocalView<Kokkos::HostSpace>();

  const size_t localNumRows = mvec->getLocalLength();
  const size_t numVectors = mvec->getNumVectors();

  Kokkos::parallel_for(localNumRows, KOKKOS_LAMBDA (const size_t row)
    {
      for (phist_lidx col = 0; col != numVectors; ++col)
      {
        reMvecView(row, col) = st::real(sourceMvecView(row, col));
        imMvecView(row, col) = st::imag(sourceMvecView(row, col));
      }
    });

  *iflag = PHIST_SUCCESS;
}
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Dconst_mvec_ptr reV, phist_Dconst_mvec_ptr imV, int *iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_MT_>::mvec_t, reMvec, reV, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_MT_>::mvec_t, imMvec, imV, *iflag);
  
  mvec->sync<Kokkos::HostSpace> ();

  auto const targetMvecView = mvec->getLocalView<Kokkos::HostSpace>();
  auto reMvecView = reMvec->getLocalView<Kokkos::HostSpace>();
  auto imMvecView = imMvec->getLocalView<Kokkos::HostSpace>();

  const size_t localNumRows = mvec->getLocalLength();
  const size_t numVectors = mvec->getNumVectors();

  Kokkos::parallel_for(localNumRows, KOKKOS_LAMBDA (const size_t row)
  {
    for (phist_lidx col = 0; col != numVectors; ++col)
    {
      targetMvecView(row, col) = reMvecView(row, col) + st::cmplx_I() * imMvecView(row, col);
    }
  });

  *iflag = PHIST_SUCCESS;
}
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Smvec* reV, phist_Smvec* imV, int *iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvec, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_MT_>::mvec_t, reMvec, reV, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_MT_>::mvec_t, imMvec, imV, *iflag);
  
  //mvec->sync<Kokkos::HostSpace> ();
  reMvec->sync<Kokkos::HostSpace> ();
  imMvec->sync<Kokkos::HostSpace> ();

  auto const sourceMvecView = mvec->getLocalView<Kokkos::HostSpace>();
  auto reMvecView = reMvec->getLocalView<Kokkos::HostSpace>();
  auto imMvecView = imMvec->getLocalView<Kokkos::HostSpace>();

  const size_t localNumRows = mvec->getLocalLength();
  const size_t numVectors = mvec->getNumVectors();

  Kokkos::parallel_for(localNumRows, KOKKOS_LAMBDA (const size_t row)
    {
      for (phist_lidx col = 0; col != numVectors; ++col)
      {
        reMvecView(row, col) = st::real(sourceMvecView(row, col));
        imMvecView(row, col) = st::imag(sourceMvecView(row, col));
      }
    });

  *iflag = PHIST_SUCCESS;
}
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Sconst_mvec_ptr reV, phist_Sconst_mvec_ptr imV, int *iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, V, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_MT_>::mvec_t, reMvec, reV, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_MT_>::mvec_t, imMvec, imV, *iflag);
  
  mvec->sync<Kokkos::HostSpace> ();

  auto const targetMvecView = mvec->getLocalView<Kokkos::HostSpace>();
  auto reMvecView = reMvec->getLocalView<Kokkos::HostSpace>();
  auto imMvecView = imMvec->getLocalView<Kokkos::HostSpace>();

  const size_t localNumRows = mvec->getLocalLength();
  const size_t numVectors = mvec->getNumVectors();

  Kokkos::parallel_for(localNumRows, KOKKOS_LAMBDA (const size_t row)
  {
    for (phist_lidx col = 0; col != numVectors; ++col)
    {
      targetMvecView(row, col) = reMvecView(row, col) + st::cmplx_I() * imMvecView(row, col);
    }
  });

  *iflag = PHIST_SUCCESS;
}
# endif
#endif

#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
//#include "../common/kernels_no_carp.cpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_io.cpp"
#include "../common/kernels_no_fused.cpp"
#include "../common/default_context_def.hpp"
